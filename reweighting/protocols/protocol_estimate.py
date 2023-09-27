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
import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params, Integer

from pwem.objects import SetOfAtomStructs, EMFile
from prody2.objects import SetOfTrajFrames
from pwchem.objects import MDSystem

import reweighting

from subprocess import check_call
import sys

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
        form.addSection(label='Inputs')

        form.addParam('infileClusterSizeData', params.EnumParam, choices=['file', 'pointer'],
                      label="Import initial cluster sizes from",
                      default=self.IMPORT_FROM_FILES,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import cluster sizes from local file or Scipion object')
        form.addParam('clusterSizeFile', params.PathParam, label="File path",
                      condition='infileClusterSizeData == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to desired cluster sizes file.')
        form.addParam('clusterSizePointer', params.PointerParam, 
                      label="Cluster size object",
                      condition='infileClusterSizeData == USE_POINTER',
                      pointerClass='SetOfAtomStructs,SetOfTrajFrames',
                      help='The input structures can be any set that contains sizes')
        
        form.addParam('infileImageDistanceData', params.EnumParam, choices=['file', 'pointer'],
                      label="Import initial cluster sizes from",
                      default=self.IMPORT_FROM_FILES,
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

        form.addParam('imageDistancePointer', params.MultiPointerParam, 
                      label="Image distance objects(s)",
                      condition='infileImageDistanceData == USE_POINTER',
                      pointerClass='EMFile',
                      help='This can be the output from a distance calculation job')
        
        form.addParam('chains', params.IntParam, default=4,
                      label="Number of MCMC chains",
                      level=params.LEVEL_ADVANCED,
                      help='Number of chains for Markov Chain Monte Carlo for posterior sampling')
        
        form.addParam('iterwarmup', params.IntParam, default=200,
                      label="Number of warm up iterations",
                      level=params.LEVEL_ADVANCED,
                      help='Number of MCMC warm up steps')
        
        form.addParam('itersample', params.IntParam, default=2000,
                      label="Number of sample iterations",
                      level=params.LEVEL_ADVANCED,
                      help='Number of MCMC sample steps')
        
        form.addParam('parallelchain', params.IntParam, default=200,
                      label="Number of chains in parallel",
                      level=params.LEVEL_ADVANCED,
                      help='(for parallelization) number of chains in parallel for MCMC')
        
        form.addParam('threadsperchain', params.IntParam, default=2000,
                      label="number of threads per chain",
                      level=params.LEVEL_ADVANCED,
                      help='(for parallelization) number of threads per chain for MCMC')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calculationStep')
        self._insertFunctionStep('createOutputStep')

    def calculationStep(self):

        if self.infileClusterSizeData.get() == self.IMPORT_FROM_FILES:
            infileclustersize = self.clusterSizeFile.get()
        else:
            infileclustersize = self.clusterSizePointer.get().getFileName()

        if self.infileImageDistanceData.get() == self.IMPORT_FROM_FILES:
            infileimagedistance = self.getMatchFiles()
        else:
            infileimagedistance = [pointer.get().getFileName() 
                                   for pointer in self.imageDistancePointer]
            
        infileimagedistance = " ".join(infileimagedistance)

        params = (infileclustersize, infileimagedistance, 
                  self._getExtraPath(),
                  self.chains.get(), self.iterwarmup.get(),
                  self.itersample.get())

        command = """python3 -m cryoER.run_cryoER_mcmc \
            --infileclustersize {0} \
            --infileimagedistance {1} --outdir {2} \
            --chains {3} \
            --iterwarmup {4} \
            --itersample {5}""".format(*params)

        command = reweighting.Plugin.getReweightingCmd(command)
        check_call(command, shell=True, stdout=sys.stdout, stderr=sys.stderr, env=None, cwd=None)

    def createOutputStep(self):
        # register output files
        self.args = {}
        filelist = sorted(os.listdir(self._getExtraPath("Stan_output")))

        for i, filename in enumerate(filelist):
            self.args["file_" + str(i+1)] = EMFile(filename=self._getExtraPath("Stan_output/" + filename))

        self._defineOutputs(**self.args)

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
