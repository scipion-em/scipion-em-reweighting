# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              James Krieger        (jamesmkrieger@gmail.com)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import numpy as np
import os

from pwem.protocols import ProtAnalysis3D
from pwem.objects import Volume, SetOfParticles

from pyworkflow import VERSION_1_1
from pyworkflow.object import Float
from pyworkflow.protocol.params import (PointerParam, StringParam, USE_GPU, GPU_LIST,
                                        BooleanParam)
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.utils as pwutils

import relion.convert

import reweighting
from reweighting.constants import REWEIGHTING_SCRIPTS

class ReweightingProtComputeLikelihood(ProtAnalysis3D):
    """This protocol computes the likelihood of a set of particles with assigned angles when compared to a
       set of maps or atomic models using CryoLike"""

    _label = 'compute likelihood'
    _lastUpdateVersion = VERSION_1_1
    _possibleOutputs = {"reprojections": SetOfParticles}

    stepsExecutionMode = STEPS_PARALLEL
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self._classesInfo = dict()

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input images", important=True,
                      pointerClass='SetOfParticles', pointerCondition='hasAlignmentProj')
        form.addParam('inputRefs', PointerParam, label="References", important=True,
                      pointerClass='Volume,SetOfVolumes,AtomStruct,SetOfAtomStructs',
                      help='Volume, set of volumes or set of atomic structures to which the set of '\
                           'particles will be compared')         
        form.addParallelSection(threads=3, mpi=8)

        form.addHidden(USE_GPU, BooleanParam, default=False,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation.\
                       Select the one you want to use. Be aware that the GPU program is new and may have problems")

        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")
    
    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()

        convId = self._insertFunctionStep(self.convertParticlesStep, prerequisites=[], needsGPU=False)
        inputRefs = self.inputRefs.get()

        stepIds = []
        i=1
        if isinstance(inputRefs, Volume):
            tempId = self._insertFunctionStep(self.convertTemplateStep, inputRefs.getFileName(), i,
                                              prerequisites=convId, needsGPU=self.useGpu)
            i += 1
            stepIds.append(tempId)
        else:
            for volume in inputRefs:
                tempId = self._insertFunctionStep(self.convertTemplateStep, volume.getFileName(), i,
                                                  prerequisites=convId, needsGPU=self.useGpu)
                i += 1
                stepIds.append(tempId)


        compStep = self._insertFunctionStep(self.compileTemplateListsStep,
                                            prerequisites=stepIds,
                                            needsGPU=False)

        stepIds = []
        i=0
        if isinstance(inputRefs, Volume):
            llId = self._insertFunctionStep(self.calculateLikelihoodStep, inputRefs.getFileName(), i,
                                            prerequisites=compStep, needsGPU=self.useGpu)
            i += 1
            stepIds.append(llId)
        else:
            for volume in inputRefs:
                llId = self._insertFunctionStep(self.calculateLikelihoodStep, volume.getFileName(), i,
                                                prerequisites=compStep, needsGPU=self.useGpu)
                i += 1
                stepIds.append(llId)

        self._insertFunctionStep(self.createOutputStep,
                                 prerequisites=stepIds,
                                 needsGPU=False)

    # --------------------------- STEPS functions ---------------------------------------------------
    def convertParticlesStep(self):
        """Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.

        Then use it to convert particles to Fourier

        Also convert image parameters like pixel size and box size to npz.
        """
        imgSet = self.inputParticles.get()
        imgStar = self._getFileName('input_star')

        outputStack = self._getExtraPath('particles/particles.mrcs')
        pwutils.makePath(self._getExtraPath("particles"))
        
        self.info("Converting set from '%s' into '%s'" %
                    (imgSet.getFileName(), imgStar))
        relion.convert.writeSetOfParticles(
            imgSet, imgStar,
            outputDir=self._getExtraPath(),
            outputStack=outputStack)

        inputParticles = self.inputParticles.get()

        args = "--input_particles_star_file %s " % imgStar
        args += '--input_particles_stack_files %s ' % outputStack
        args += '--folder_output %s ' % self._getExtraPath()

        Ts = inputParticles.getSamplingRate()
        xdim = inputParticles.getFirstItem().getXDim()
        args += '--pixel_size %f --box_size %d ' % (Ts, xdim)

        args += '--batch_size %d ' % inputParticles.getSize()
        if self.useGpu:
            args+="--use_cuda "
            gpuId = self._stepsExecutor.getGpuList()
            if isinstance(gpuId, int):
                gpuStr = str(gpuId)
            else:
                gpuStr = ','.join([str(g) for g in gpuId])
            os.environ["CUDA_VISIBLE_DEVICES"] = gpuStr

        command = "python3 " + os.path.join(REWEIGHTING_SCRIPTS, "cryolike_inputs.py")
        self.runJob(reweighting.Plugin.getCryoLikeCmd(command), args)

    def convertTemplateStep(self, fnVol, i):
        """Create Fourier reprojection templates for this volume and use converted them
        with converted particles and image parameters to calculate likelihoods.
        """
        args = '--i %d --folder_output %s --ref %s ' % (i, self._getExtraPath(), fnVol)
        if self.useGpu:
            args+="--use_cuda "
            gpuId = self._stepsExecutor.getGpuList()
            if isinstance(gpuId, int):
                gpuStr = str(gpuId)
            else:
                gpuStr = ','.join([str(g) for g in gpuId])
            os.environ["CUDA_VISIBLE_DEVICES"] = gpuStr
        command = "python3 " + os.path.join(REWEIGHTING_SCRIPTS, "cryolike_template.py")
        self.runJob(reweighting.Plugin.getCryoLikeCmd(command), args)


    def compileTemplateListsStep(self):
        """Read template list npy files and compile them into one list npy file.
        """
        list_of_file_lists = pwutils.glob(self._getExtraPath('templates/*npy'))
        list_of_template_files = []
        for list_file in list_of_file_lists:
            list_of_template_files.extend(list(np.load(list_file)))

        np.save(self._getExtraPath('templates/template_file_list.npy'), 
                list_of_template_files)

    def calculateLikelihoodStep(self, fnVol, i):
        """Create Fourier reprojection templates for this volume and use converted them
        with converted particles and image parameters to calculate likelihoods.
        """
        args = '--i %d --folder_output %s --ref %s ' % (i, self._getExtraPath(), fnVol)
        if self.useGpu:
            args+="--use_cuda "
            gpuId = self._stepsExecutor.getGpuList()
            if isinstance(gpuId, int):
                gpuStr = str(gpuId)
            else:
                gpuStr = ','.join([str(g) for g in gpuId])
            os.environ["CUDA_VISIBLE_DEVICES"] = gpuStr
        command = "python3 " + os.path.join(REWEIGHTING_SCRIPTS, "cryolike_likelihood.py")
        self.runJob(reweighting.Plugin.getCryoLikeCmd(command), args)        

    def createOutputStep(self):
        inputParticles = self.inputParticles.get()

        outputSet = self._createSetOfParticles()
        outputSet.copyInfo(inputParticles)

        self.LLs = {}

        self.idmap = {}
        for j, ind in enumerate([item.getObjId() for item in inputParticles]):
            self.idmap[ind] = j

        self.i=1
        if isinstance(self.inputRefs.get(), Volume):
            self.LLs[self.i] = np.load(os.path.join(self._getExtraPath(
                'likelihood/template%d/log_likelihood/log_likelihood_integrated_fourier_stack_000000.npy') % self.i))
            outputSet.copyItems(self.inputParticles.get(), updateItemCallback=self._processRow)
            self.i += 1
        else:
            for _ in self.inputRefs.get():
                self.LLs[self.i] = np.load(os.path.join(self._getExtraPath(
                    'likelihood/template%d/log_likelihood/log_likelihood_integrated_fourier_stack_000000.npy') % self.i))
                outputSet.copyItems(self.inputParticles.get(), updateItemCallback=self._processRow)
                self.i += 1

        self._defineOutputs(reprojections=outputSet)
        self._defineSourceRelation(self.inputParticles, outputSet)

        matrix = np.array([particle._cryolike_logLikelihood.get() for particle in outputSet])
        matrix = matrix.reshape((self.i-1,-1))
        np.save(self._getExtraPath('matrix.npy'), matrix)

    def _processRow(self, particle, row):
        setattr(particle, '_cryolike_logLikelihood', 
                Float(self.LLs[self.i][self.idmap[particle.getObjId()]]))
        particle.setObjId(None)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_star': self._getPath('input_particles.star'),
        }
        self._updateFilenamesDict(myDict)
