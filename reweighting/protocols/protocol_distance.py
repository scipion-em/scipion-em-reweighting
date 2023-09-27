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
This module will provide the image distances program
"""
import reweighting

import os

from pwem.objects import AtomStruct, SetOfTrajFrames, SetOfParticles, EMFile
from pwem.protocols import EMProtocol

from pwchem.objects import MDSystem

from pyworkflow.protocol import params, Integer

from subprocess import check_call
import sys

class ReweightingImageDistancesProtocol(EMProtocol):
    """
    This protocol will calculate image distances from structures and particles
    """
    _label = 'Calc image distances'

    IMPORT_FROM_FILES = 0
    USE_POINTER = 1

    USE_ATOMSTRUCT = 1
    USE_MDSYSTEM = 2

    CUDA = 0
    CPU = 1

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label='Input structures')

        form.addParam('topImageDataType', params.EnumParam, choices=['file', 'atomstruct', 'mdsystem'],
                      label="Import images topology from",
                      default=self.USE_ATOMSTRUCT,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Choose whether to import images topology from local file or using a pointer')
        form.addParam('topImageFile', params.PathParam, label="Images topology file path",
                      condition='topImageDataType == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to images topology.')
        form.addParam('topImageStructure', params.PointerParam, label="Images topology/structure",
                      condition='topImageDataType == USE_ATOMSTRUCT or topImageDataType == USE_MDSYSTEM',
                      pointerClass='AtomStruct, MDSystem',
                      help="The images topology can currently be provided as an AtomStruct or MDSystem. "
                           "If it's an MDSystem then it also provides the trajectory.")
        
        form.addParam('trajImageDataType', params.EnumParam, choices=['file', 'pointer'],
                      label="Import images trajectory from",
                      default=self.USE_POINTER,
                      condition='topImageDataType != USE_MDSYSTEM',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import images trajectory from local file or using a pointer')
        form.addParam('trajImageFile', params.PathParam, label="Images trajectory file path",
                      condition='trajImageDataType == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to images trajectory.')
        form.addParam('trajImage', params.PointerParam, label="Images trajectory",
                      condition='trajImageDataType == USE_ATOMSTRUCT and topImageDataType != USE_MDSYSTEM',
                      pointerClass='AtomStruct, SetOfTrajFrames',
                      help='The images trajectory can currently be provided as an AtomStruct, MDSystem or SetOfTrajFrames')
        
        form.addParam('topStructDataType', params.EnumParam, choices=['file', 'atomstruct', 'mdsystem'],
                      label="Import topology atomic structure from",
                      default=self.USE_POINTER,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Whether to import structures topology from local file or using a pointer')
        form.addParam('topStructFile', params.PathParam, label="Structures topology file path",
                      condition='topStructDataType == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to topology struc.')
        form.addParam('topStructureStructure', params.PointerParam, label="Structures topology/structure",
                      condition='topStructDataType == USE_ATOMSTRUCT or topStructDataType == USE_MDSYSTEM',
                      pointerClass='AtomStruct, MDSystem',
                      help="The structures topology can currently be provided as an AtomStruct or MDSystem. "
                           "If it's an MDSystem then it also provides the trajectory.")
        
        form.addParam('trajStructDataType', params.EnumParam, choices=['file', 'pointer'],
                      label="Import structures trajectory from",
                      default=self.USE_POINTER,
                      condition='topStructDataType != USE_MDSYSTEM',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import structures trajectory from local file or using a pointer')
        form.addParam('trajStructFile', params.PathParam, label="Structures trajectory file path",
                      condition='trajStructDataType == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to images trajectory.')
        form.addParam('trajStruct', params.PointerParam, label="Structures trajectory",
                      condition='trajStructDataType == USE_ATOMSTRUCT and topStructDataType != USE_MDSYSTEM',
                      pointerClass='AtomStruct, SetOfTrajFrames',
                      help='The structures trajectory can currently be provided as an AtomStruct, MDSystem or SetOfTrajFrames')

        form.addParam('device', params.EnumParam, choices=['cuda', 'cpu'],
                      label="hardware device",
                      default=self.CPU,
                      expertLevel=params.LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='hardware device for calculation: "cuda" for GPU, or "cpu" for CPU')

        form.addParam('nBatch', params.IntParam, default=10,
                      label="Number of batches",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='number of batches to separate the output files into '
                           'for memory management')


        form.addSection(label='Simulated image parameters')

        form.addParam('nPixel', params.IntParam, default=128,
                      label="Number of pixels",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='number of image pixels, use power of 2 for CTF purpose')

        form.addParam('pixelSize', params.FloatParam, default=1.,
                      label="Pixel size in Angstrom",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Pixel size in Angstrom')
        
        form.addParam('sigma', params.FloatParam, default=1.5,
                      label="radius of Gaussian atom",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='radius of Gaussian for atoms')
        
        form.addParam('snr', params.FloatParam, default=1e-2,
                      label="signal-to-noise ratio",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='signal-to-noise ratio')
        
        form.addParam('ctfBool', params.BooleanParam, default=False,
                      label="introduce CTF modulation?",
                      expertLevel=params.LEVEL_ADVANCED,
                      help='whether to introduce CTF modulation')
        

        
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calcDistanceStep')
        self._insertFunctionStep('createOutputStep')

    def calcDistanceStep(self):
        if self.topImageDataType.get() == self.IMPORT_FROM_FILES:
            topImageFile = self.topImageFile.get()
        elif self.topImageDataType.get() == self.USE_ATOMSTRUCT:
            topImageFile = self.topImageStructure.get().getFileName()
        else:
            imageSystem = self.topImageStructure.get()
            topImageFile = imageSystem.getTopologyFile()

        if self.trajImageDataType.get() == self.IMPORT_FROM_FILES:
            trajImageFile = self.trajImageFile.get()
        elif (self.trajImageDataType.get() == self.USE_ATOMSTRUCT 
              and self.topImageDataType != self.USE_MDSYSTEM):
            trajImageFile = self.trajImage.get().getFileName()
        else:
            trajImageFile = imageSystem.getTrajectoryFile()

        if self.topStructDataType.get() == self.IMPORT_FROM_FILES:
            topStructFile = self.topStructFile.get()
        elif self.topStructDataType.get() == self.USE_ATOMSTRUCT:
            topStructFile = self.topStructureStructure.get().getFileName()
        else:
            structSystem = self.topImageStructure.get()
            topStructFile = structSystem.getTopologyFile()        

        if self.trajStructDataType.get() == self.IMPORT_FROM_FILES:
            trajStructFile = self.trajStructFile.get()
        elif (self.trajStructDataType.get() == self.USE_ATOMSTRUCT 
              and self.topStructDataType != self.USE_MDSYSTEM):
            trajStructFile = self.trajStructureStructure.get().getFileName()
        else:
            trajStructFile = structSystem.getTrajectoryFile()            

        matricesFilename = self._getExtraPath('rot_mats_struc_image.npy')

        if self.device.get() == self.CPU:
            device = "cpu"
        else:
            device = "cuda"

        params = (topImageFile, trajImageFile, 
                  topStructFile, trajStructFile,
                  matricesFilename, self._getExtraPath(),
                  self.nPixel.get(), self.pixelSize.get(),
                  self.sigma.get(), self.snr.get(),
                  self.nBatch.get(), device)
        
        command = """python3 -m tools.calc_rot_mats --top_image {0} --traj_image {1} \
            --top_struc {2} --traj_struc {3} --outdir {5} --n_batch 5""".format(*params)
        command = reweighting.Plugin.getReweightingCmd(command)
        check_call(command, shell=True, stdout=sys.stdout, stderr=sys.stderr, env=None, cwd=None)

        command = """python3 -m cryoER.calc_image_struc_distance --top_image {0} --traj_image {1} \
            --top_struc {2} --traj_struc {3} --rotmat_struc_imgstruc {4} \
            --outdir {5} \
            --n_pixel {6} \
            --pixel_size {7} \
            --sigma {8} \
            --signal_to_noise_ratio {9} \
            --n_batch {10} --device {11}""".format(*params)
        
        if self.ctfBool.get():
            command += " --ctf"

        command = reweighting.Plugin.getReweightingCmd(command)
        check_call(command, shell=True, stdout=sys.stdout, stderr=sys.stderr, env=None, cwd=None)

    def createOutputStep(self):
        # register output files
        self.args = {}
        for i, filename in enumerate(os.listdir(self._getExtraPath())):
            self.args["File" + str(i+1)] = EMFile(filename=filename)

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
            methods.append("Protocol has been printed methods.")
        return methods
