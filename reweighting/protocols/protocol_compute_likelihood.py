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
try:
    from itertools import izip
except ImportError:
    izip = zip

import numpy as np
import os

from pwem.protocols import ProtAnalysis3D
from pwem.objects import Volume, SetOfParticles, SetOfClasses3D

from pyworkflow import VERSION_1_1
from pyworkflow.object import Float
from pyworkflow.protocol.params import (PointerParam, StringParam, USE_GPU, GPU_LIST,
                                        BooleanParam, FloatParam, IntParam)
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.utils as pwutils

import relion.convert

import reweighting
from reweighting.constants import REWEIGHTING_SCRIPTS

class ReweightingProtComputeLikelihood(ProtAnalysis3D):
    """
    Computes the likelihood of a set of particles with assigned projection
    angles against one or more reference maps or atomic models using
    CryoLike.

    AI Generated:

    Compute Likelihood (ReweightingProtComputeLikelihood) — User Manual
        Overview

        This protocol evaluates how well a set of experimental particle
        images agrees with one or more structural references by computing
        likelihood scores in Fourier space. Its main objective is to
        quantify, for each particle, the statistical compatibility between
        the observed image and a collection of reference projections
        generated from maps or atomic structures.

        In cryo-EM workflows, this type of likelihood estimation is
        especially useful when particles already carry assigned projection
        orientations and the user wishes to compare how strongly different
        structural hypotheses explain the experimental data. Typical
        applications include reweighting heterogeneous datasets, assigning
        particles to competing structural states, or evaluating how well
        candidate references capture conformational variability.

        Inputs and General Workflow

        The protocol requires a set of aligned particles and one or more
        structural references. The particle images must already contain
        projection alignment information, since the likelihood calculation
        relies on known viewing geometry.

        During execution, the protocol first converts the input particles
        into the internal format expected by CryoLike. Particle metadata,
        Fourier-transformed particle images, and acquisition parameters
        such as pixel size and box dimensions are prepared for downstream
        computation.

        Each reference is then processed independently. For every map or
        atomic model, CryoLike generates reprojection templates in Fourier
        space covering a sampled set of orientations and in-plane
        rotations. These templates form the reference ensemble against
        which the experimental particles are evaluated.

        Sampling of Orientation and Image Displacements

        The protocol samples projection space using two complementary
        components.

        The first component is angular sampling. A viewing distance
        controls how densely orientations are distributed on the sphere,
        while the number of in-plane rotations defines rotational sampling
        around each projection direction.

        The second component is translational sampling. The protocol
        explores possible image shifts by defining a maximum displacement
        in pixels and subdividing that range into a grid of X and Y
        displacements.

        From a practical cryo-EM perspective, broader sampling increases
        robustness when alignment uncertainty is high, whereas narrower
        sampling can significantly reduce computational cost when
        orientations are already well refined.

        Likelihood Computation

        Once particles and templates are prepared, CryoLike computes
        integrated log-likelihood values for every particle-reference
        combination.

        These values represent how well each reference explains each
        experimental particle after integrating over sampled orientations,
        in-plane rotations, and translational displacements. The resulting
        likelihoods provide a quantitative basis for comparing multiple
        structural hypotheses.

        When several references are supplied, the protocol evaluates all
        of them independently and stores the likelihood values in a
        particle-by-reference matrix. This matrix can later be used for
        downstream statistical analysis, classification refinement, or
        likelihood-based particle reweighting.

        CPU and GPU Execution

        The protocol supports both CPU and GPU execution.

        On CPU systems, parallel execution is controlled by Scipion-level
        parallelization together with the number of threads used internally
        by CryoLike. Since several CryoLike calls may run simultaneously,
        the effective resource usage depends on both levels of
        parallelization.

        On GPU systems, the protocol can distribute CryoLike execution
        across selected CUDA devices. This can substantially accelerate
        template generation and likelihood computation, especially when
        processing large particle datasets or multiple references.

        In practice, memory availability becomes an important limiting
        factor because particle batches, Fourier templates, and likelihood
        arrays may all need to reside in memory at the same time.

        Outputs and Interpretation

        The protocol produces an output particle set in which every
        particle carries a CryoLike log-likelihood value.

        When multiple references are provided, the protocol also compiles
        a likelihood matrix whose rows correspond to references and whose
        columns correspond to particles. This matrix provides a compact
        quantitative summary of the agreement between experimental data
        and structural hypotheses.

        Based on the maximum likelihood value for each particle, the
        protocol automatically creates a 3D classification output. Each
        particle is assigned to the reference that best explains it.

        Biologically, this output can be interpreted as a likelihood-based
        partitioning of particles among alternative structural states,
        making it particularly useful for heterogeneous systems,
        conformational landscapes, or competing atomic models.

        Practical Considerations

        The biological meaning of the likelihood values depends strongly
        on the quality and relevance of the references. Closely related
        references may produce subtle likelihood differences, whereas
        highly distinct structures often lead to clearer separation.

        The quality of the assigned particle orientations also strongly
        influences reliability. If projection parameters are poor, the
        likelihood scores may become less discriminative or biologically
        ambiguous.

        For exploratory analyses, moderate angular and translational
        sampling often provides a good balance between runtime and
        robustness. For publication-level analyses or difficult
        heterogeneous datasets, denser sampling may improve reliability at
        the cost of substantially increased computation.

        Final Perspective

        This protocol turns a set of already aligned particle images into
        a quantitative comparison against one or more structural
        references.

        Rather than performing reconstruction directly, it measures how
        strongly each particle supports competing structural models. In
        cryo-EM reweighting and heterogeneity analysis, this provides a
        statistically meaningful bridge between experimental images and
        structural interpretation.
    """

    _label = 'compute likelihood'
    _lastUpdateVersion = VERSION_1_1
    _possibleOutputs = {"reprojections": SetOfParticles}

    stepsExecutionMode = STEPS_PARALLEL
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self._classesInfo = dict()

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addParam('binThreads', IntParam,
                      label='threads',
                      default=2,
                      help='Number of threads used by CryoLike each time it is called in the protocol execution. For '
                           'example, if 3 Scipion threads and 3 CryoLike threads are set, the particles will be '
                           'processed in groups of 2 at the same time with a call of CryoLike with 3 threads each, so '
                           '6 threads will be used at the same time. Beware the memory of your machine has '
                           'memory enough to load together the number of particles specified by Scipion threads.')

        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input images", important=True,
                      pointerClass='SetOfParticles', pointerCondition='hasAlignmentProj')
        form.addParam('inputRefs', PointerParam, label="References", important=True,
                      pointerClass='Volume,SetOfVolumes,AtomStruct,SetOfAtomStructs',
                      help='Volume, set of volumes or set of atomic structures to which the set of '\
                           'particles will be compared')

        form.addParam('viewingDistance', FloatParam, label="Viewing Distance", default=1,
                      help='Distance between views sampled in units of 1/4*pi on the sphere of rotations, '\
                           'excluding in-plane ones')
        form.addParam('nInplanes', IntParam, label="Number of In-Plane Rotations", default=256,
                      help='Number of in-plane rotations sampled for each view')

        form.addParam('max_displacement_pixels', IntParam, label="Maximum displacement in pixels", default=8,
                      help='Distance between views sampled in units of 1/4*pi on the sphere of rotations, '\
                           'excluding in-plane ones')
        form.addParam('n_displacements_x', IntParam, label="Number of displacements in x", default=16,
                      help='The maximum displacement is divided into this many displacements')
        form.addParam('n_displacements_y', IntParam, label="Number of displacements in y", default=16,
                      help='The maximum displacement is divided into this many displacements')

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

        args += '--viewing_distance %f --n_inplanes %d ' % (self.viewingDistance.get(),
                                                            self.nInplanes.get())

        args += '--batch_size %d ' % len(inputParticles)
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
        list_of_file_lists = sorted(pwutils.glob(self._getExtraPath('templates/*npy')))
        list_of_template_files = []
        for list_file in list_of_file_lists:
            filename = os.path.split(list_file)[1]
            if (filename[:-4].replace('_','').isalnum()
                and not filename[:-4].replace('_','').isalpha()):
                list_of_template_files.extend(list(np.load(list_file)))

        np.save(self._getExtraPath('templates/template_file_list.npy'),
                list_of_template_files)

    def calculateLikelihoodStep(self, fnVol, i):
        """Create Fourier reprojection templates for this volume and use converted them
        with converted particles and image parameters to calculate likelihoods.
        """
        args = '--i %d --folder_output %s --ref %s ' % (i, self._getExtraPath(), fnVol)
        args += '--batch_size %d --template_batch_size %d ' % (len(self.inputParticles.get()), 8)
        args += '--max_displacement_pixels %d --n_displacements_x %d --n_displacements_y %d ' % (self.max_displacement_pixels.get(),
                                                                                                 self.n_displacements_x.get(),
                                                                                                 self.n_displacements_y.get())
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

        refsDict = {}
        self.i=1
        if isinstance(self.inputRefs.get(), Volume):
            self.LLs[self.i] = np.load(os.path.join(self._getExtraPath(
                'likelihood/template%d/log_likelihood/log_likelihood_integrated_fourier_stack_000000.npy') % (self.i-1)))
            outputSet.copyItems(self.inputParticles.get(), updateItemCallback=self._processRow)
            refsDict[self.i] = self.inputRefs.get()
            self.i += 1
        else:
            for item in self.inputRefs.get():
                self.LLs[self.i] = np.load(os.path.join(self._getExtraPath(
                    'likelihood/template%d/log_likelihood/log_likelihood_integrated_fourier_stack_000000.npy') % (self.i-1)))
                outputSet.copyItems(self.inputParticles.get(), updateItemCallback=self._processRow)
                refsDict[self.i] = item.clone()
                self.i += 1

        self._defineOutputs(reprojections=outputSet)
        self._defineSourceRelation(self.inputParticles, outputSet)

        matrix = np.array([particle._cryolike_logLikelihood.get() for particle in outputSet])
        matrix = matrix.reshape((self.i-1,-1))
        np.save(self._getExtraPath('matrix.npy'), matrix)

        classIds = np.argmax(matrix, axis=0)+1

        clsSet = SetOfClasses3D.create(self._getExtraPath())
        clsSet.setImages(inputParticles)

        clsDict = {}  # Dictionary to store the (classId, classSet) pairs

        for ref, rep in refsDict.items():
            # add empty classes
            classItem = clsSet.ITEM_TYPE.create(self._getExtraPath(), suffix=ref+1)
            classItem.setRepresentative(rep)
            clsDict[ref] = classItem
            clsSet.append(classItem)

        cls_prev = 1
        for img, ref in izip(inputParticles, classIds):
            if ref != cls_prev:
                cls_prev = ref

            classItem = clsDict[ref]
            classItem.append(img)

        for classItem in clsDict.values():
            clsSet.update(classItem)

        clsSet.write()

        self._defineOutputs(outputClasses=clsSet)
        self._defineSourceRelation(self.inputParticles, clsSet)
        self._defineSourceRelation(self.inputRefs, clsSet)

    def _processRow(self, particle, row):
        setattr(particle, '_cryolike_logLikelihood', 
                Float(self.LLs[self.i][self.idmap[particle.getObjId()]]))
        particle.setObjId(None)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_star': self._getExtraPath('particles/input_particles.star'),
        }
        self._updateFilenamesDict(myDict)
