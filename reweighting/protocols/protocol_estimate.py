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
import numpy as np
import os

from pwem.protocols import EMProtocol
from pwem.objects import EMSet, EMObject, EMFile

from pyworkflow.protocol import params
from pyworkflow.object import Float, String

import reweighting
from reweighting.constants import (REWEIGHTING_SCRIPTS, 
                                   REWEIGHTING_MEAN, REWEIGHTING_STD)

class ReweightingEstimateWeightsProtocol(EMProtocol):
    """
    This protocol will estimate weights using MCMC
    """
    _label = 'Estimate weights'

    IMPORT_FROM_FILES = 0
    USE_POINTER = 1

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addParallelSection(threads=1, mpi=0)

        form.addSection(label='Inputs')

        form.addParam('infileClusterSizeData', params.EnumParam, choices=['file', 'pointer'],
                      label="Import initial cluster sizes from",
                      default=self.USE_POINTER,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import cluster sizes from local file or Scipion object')
        form.addParam('clusterSizeFile', params.PathParam, label="File path",
                      condition='infileClusterSizeData == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to desired cluster sizes file.')
        form.addParam('clusterSizePointer', params.PointerParam, 
                      label="Cluster size object",
                      condition='infileClusterSizeData == USE_POINTER',
                      pointerClass='SetOfAtomStructs,SetOfTrajFrames,SetOfVolumes,SetOfClassesTraj',
                      help='The input structures can be any set that contains sizes')
        
        form.addParam('infileImageDistanceData', params.EnumParam, choices=['file', 'pointer'],
                      label="Import log likelihood or distances from",
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
                      pointerClass='EMFile,SetOfParticles',
                      help='This can be the output from a distance calculation job')
        
        form.addParam('chains', params.IntParam, default=4,
                      label="Number of MCMC chains",
                      help='Number of chains for Markov Chain Monte Carlo for posterior sampling')
        
        form.addParam('iterwarmup', params.IntParam, default=200,
                      label="Number of warm up iterations",
                      help='Number of MCMC warm up steps')
        
        form.addParam('itersample', params.IntParam, default=2000,
                      label="Number of sample iterations",
                      help='Number of MCMC sample steps')
        
        form.addParam('parallelchain', params.IntParam, default=1,
                      label="Number of chains in parallel",
                      help='(for parallelization) number of chains in parallel for MCMC')

        form.addParam('lambda_', params.FloatParam, default=-1,
                      label="Noise standard deviation lambda",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Leave it as -1 if using logLikelihood to not normalise with lambda again')  

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('calculationStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):

        if self.infileClusterSizeData.get() == self.IMPORT_FROM_FILES:
            self.infileclustersize = self.clusterSizeFile.get()
        else:
            clusterSizeObject = self.clusterSizePointer.get()
            if hasattr(clusterSizeObject.getFirstItem(), '_weights'):
                self.clusterSizes = [item._weights.get() for item in clusterSizeObject]
            elif hasattr(clusterSizeObject.getFirstItem(), '_prodyWeights'):
                self.clusterSizes = [item._prodyWeights.get() for item in clusterSizeObject]
            else:
                self.clusterSizes = np.ones(len(clusterSizeObject))

            if isinstance(self.clusterSizes, list):
                self.clusterSizes = np.array(self.clusterSizes)

            if self.clusterSizes.sum() != 1:
                self.clusterSizes /= self.clusterSizes.sum()

            self.infileclustersize = self._getExtraPath('cluster_sizes.txt')
            np.savetxt(self.infileclustersize, self.clusterSizes)

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

                infileimagedistance.append(filename)
            
        self.infileimagedistance = " ".join(infileimagedistance)

    def calculationStep(self):

        params = (self.infileclustersize, self.infileimagedistance, 
                  self._getExtraPath(),
                  self.chains.get(), self.iterwarmup.get(),
                  self.itersample.get(), self.lambda_.get())

        command = "python3 -m cryoER.run_cryoER_mcmc"
        args = """--infileclustersize {0} --infileimagedistance {1} --outdir \
{2} --chains {3} --iterwarmup {4} --itersample {5} --lmbd {6}""".format(*params)
        
        parallelChains = self.parallelchain.get()
        threadsPerChain = int(self.numberOfThreads.get()/parallelChains)
        if parallelChains > 1 or threadsPerChain > 1:
            args += " --parallelchain {0} --threadsperchain {1}".format(parallelChains, 
                                                                        threadsPerChain)

        self.runJob(reweighting.Plugin.getReweightingCmd(command), args)

        command2 = "python3 " + os.path.join(REWEIGHTING_SCRIPTS, "analyse.py")
        args2 = "--output_directory {0} --filename_cluster_counts {1}".format(self._getExtraPath(),
                                                                              self.infileclustersize)
        self.runJob(reweighting.Plugin.getReweightingCmd(command2), args2)

    def createOutputStep(self):
        # register output files
        self.args = {}

        self.means = np.loadtxt(self._getExtraPath("reweighting_mean_weights.txt"))
        self.stds = np.loadtxt(self._getExtraPath("reweighting_std_weights.txt"))

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

        self.args["outputSet"] = outSet
        self._defineOutputs(**self.args)

    def _addWeights(self, item, row=None):
        idx = self.idxMap.index(item.getObjId())

        # We provide data directly so don't need a row
        mean = Float(self.means[idx])
        mean.setPrecision(1e-6)
        setattr(item, REWEIGHTING_MEAN, mean)

        std = Float(self.stds[idx])
        std.setPrecision(1e-6)
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
