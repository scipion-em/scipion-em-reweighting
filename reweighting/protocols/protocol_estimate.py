# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *
# * Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


"""
This module will estimate weights using MCMC
"""
from multiprocessing import cpu_count
import numpy as np
import os

from pwem.objects import EMFile, EMSet, EMObject
from pwem.protocols import EMProtocol

from pyworkflow.protocol import params
from pyworkflow.object import Float

import reweighting
from reweighting.constants import (REWEIGHTING_SCRIPTS, 
                                   REWEIGHTING_MEAN, REWEIGHTING_STD)

class ReweightingEstimateWeightsProtocol(EMProtocol):
    """
    This protocol will estimate weights using MCMC
    """
    _label = 'Reweight'

    USE_POINTER = 0
    USE_TEXT = 1
    IMPORT_FROM_FILES = 1

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=0)

        form.addSection(label='Inputs')

        form.addParam('infileClusterSizeData', params.EnumParam, choices=['pointer', 'text'],
                      label="Import initial cluster sizes from",
                      default=self.USE_POINTER,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Read cluster sizes from Scipion Set object or enter your own text')
        form.addParam('clusterSizePointer', params.PointerParam, 
                      label="Cluster size object",
                      condition='infileClusterSizeData == USE_POINTER',
                      pointerClass='SetOfAtomStructs,SetOfTrajFrames',
                      help='The input structures can be any set that contains weights. '
                           'If there are none, then the prior is assumed to be uniform.')
        form.addParam('self.clusterSizes', params.StringParam,
                      label="Cluster sizes",
                      condition='infileClusterSizeData == USE_TEXT',
                      help='Enter prior cluster sizes.')
        
        form.addParam('infileImageDistanceData', params.EnumParam, choices=['pointer', 'file'],
                      label="Import maximum likelihood distances from",
                      default=self.USE_POINTER,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import cluster sizes from local file or Scipion object')

        form.addParam('filesPath', params.PathParam,
                      condition="infileImageDistanceData == IMPORT_FROM_FILES",
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/project/data/day??_files/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/project/data/day*_files/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")

        form.addParam('filesPattern', params.StringParam,
                      label='File pattern',
                      condition="infileImageDistanceData == IMPORT_FROM_FILES",
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters" 
                           "('*', '?', '#', ':', '%') cannot appear in the actual path.\n\n")

        form.addParam('imageDistancePointers', params.MultiPointerParam, 
                      label="Image distance object(s)",
                      condition='infileImageDistanceData == USE_POINTER',
                      pointerClass='EMFile, SetOfParticles',
                      help='This can be the output from a distance or logLikelihood calculation job')
        
        form.addParam('lambda_', params.FloatParam, default=np.sqrt(0.5),
                      label="Noise standard deviation lambda",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Leave it as -sqrt(0.5) if using logLikelihood to not normalise')        
        
        form.addParam('chains', params.IntParam, default=4,
                      label="Number of MCMC chains",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Number of chains for Markov Chain Monte Carlo for posterior sampling')
        
        form.addParam('iterwarmup', params.IntParam, default=200,
                      label="Number of warm up iterations",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Number of MCMC warm up steps')
        
        form.addParam('itersample', params.IntParam, default=2000,
                      label="Number of sample iterations",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Number of MCMC sample steps')
        
        form.addParam('parallelchain', params.IntParam, default=1,
                      label="Number of chains in parallel",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='(for parallelization) number of chains in parallel for MCMC')
        
        form.addParam('threadsperchain', params.IntParam, default=1,
                      label="number of threads per chain",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='(for parallelization) number of threads per chain for MCMC')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calculationStep')
        self._insertFunctionStep('createOutputStep')

    def calculationStep(self):

        if self.infileClusterSizeData.get() == self.USE_POINTER:
            clusterSizeObject = self.clusterSizePointer.get()
            if hasattr(clusterSizeObject.getFirstItem(), '_weights'):
                self.clusterSizes = [item._weights.get() for item in clusterSizeObject]
            elif hasattr(clusterSizeObject.getFirstItem(), '_prodyWeights'):
                self.clusterSizes = [item._prodyWeights.get() for item in clusterSizeObject]
            else:
                self.clusterSizes = np.ones(len(clusterSizeObject))
        else:
            self.clusterSizes = self.self.clusterSizes.get()

        if isinstance(self.clusterSizes, list):
            self.clusterSizes = np.array(self.clusterSizes)

        if self.clusterSizes.sum() != 1:
            self.clusterSizes /= self.clusterSizes.sum()

        infileclustersize = self._getExtraPath('cluster_sizes.txt')
        np.savetxt(infileclustersize, self.clusterSizes)

        # logger.info(redStr("prior cluster sizes {0} written to {1}".format(self.clusterSizes, infileclustersize)))

        if self.infileImageDistanceData.get() == self.IMPORT_FROM_FILES:
            infileimagedistance = self.getMatchFiles()
        else:
            infileimagedistance = []
            for i, pointer in enumerate(self.imageDistancePointers):
                distanceObject = pointer.get()
                if isinstance(distanceObject, EMFile):
                    filename = distanceObject.getFileName()
                else:
                    if hasattr(distanceObject[1], '_xmipp_logLikelihood'):
                        imageDistances = np.array([particle._xmipp_logLikelihood.get() for particle in distanceObject])
                        imageDistances = imageDistances.reshape((len(self.clusterSizes),-1))
                        filename = self._getExtraPath('image_distances_{0}.npy'.format(i+1))
                        np.save(filename, imageDistances)
                        # logger.info(redStr("image distances with shape {0} written to {1}".format(imageDistances.shape, filename)))
                
                infileimagedistance.append(filename)
        
        infileimagedistance = " ".join(infileimagedistance)

        params = (infileclustersize, infileimagedistance, 
                  self._getExtraPath(),
                  self.chains.get(), self.iterwarmup.get(),
                  self.itersample.get(), self.lambda_.get())

        command = "python3 -m cryoER.run_cryoER_mcmc"
        args = """--infileclustersize {0} --infileimagedistance {1} --outdir {2} \
--chains {3} --iterwarmup {4} --itersample {5} --lmbd {6}""".format(*params)
        
        parallelChains = self.parallelchain.get()
        threadsPerChain = self.threadsperchain.get()
        if parallelChains > 1 or threadsPerChain > 1:
            args += " --parallelchain {0} --threadsperchain {1}".format(parallelChains, 
                                                                        threadsPerChain)

        self.runJob(reweighting.Plugin.getReweightingCmd(command), args)

        command2 = "python3 " + os.path.join(REWEIGHTING_SCRIPTS, "analyse.py")
        args2 = "--output_directory {0} --filename_cluster_counts {1}".format(self._getExtraPath(),
                                                                              infileclustersize)
        self.runJob(reweighting.Plugin.getReweightingCmd(command2), args2)

    def createOutputStep(self):
        self.means = np.loadtxt(self._getExtraPath("mean_weights.txt"))
        self.stds = np.loadtxt(self._getExtraPath("std_weights.txt"))
        
        args = {}
        inSet = self.clusterSizePointer.get()
        if inSet is not None:
            inputClass = type(inSet)
            self.idxMap = list(inSet.getIdSet())
            outSet = inputClass().create(self._getExtraPath())
            outSet.copyItems(inSet, updateItemCallback=self._addWeights)
        else:
            outSet = EMSet().create(self._getExtraPath())
            self.idxMap = []
            for i, _ in enumerate(self.means):
                self.idxMap.append(i)
                item = EMObject()
                self._addWeights(item)
                outSet.append(item)

        args["outputSet"] = outSet

        self._defineOutputs(**args)

    def _addWeights(self, item, row=None):
        idx = self.idxMap.index(item.getObjId())

        # We provide data directly so don't need a row
        mean = Float(self.means[idx])
        setattr(item, REWEIGHTING_MEAN, mean)

        std = Float(self.stds[idx])
        setattr(item, REWEIGHTING_STD, std)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        if self.isFinished():
            summary.append("This protocol has finished.")
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            methods.append("This protocol has printed methods")
        return methods
