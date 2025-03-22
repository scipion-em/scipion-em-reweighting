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
This module will calculate correlations between log likelihood matrices
"""
import numpy as np
import os

from pwem.protocols import EMProtocol
from pwem.objects import EMSet, EMObject, EMFile

from pyworkflow.protocol.params import PointerParam, BooleanParam
from pyworkflow.object import Float

class ReweightingCorrelateProtocol(EMProtocol):
    """
    This protocol will calculate correlations between log likelihood matrices
    """
    _label = 'Correlate log likelihoods'

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

        form.addParam('inputParticles1', PointerParam, label="Input images 1", important=True,
                      pointerClass='SetOfParticles', pointerCondition='hasAlignmentProj')

        form.addParam('inputParticles2', PointerParam, label="Input images 2", important=True,
                      pointerClass='SetOfParticles', pointerCondition='hasAlignmentProj')

        form.addParam('flipVols', BooleanParam, label="Flip volumes for matrix 2?", default=False,
                      help='Select whether to flip volumes for matrix 2. This may be useful '
                           'if they are in the wrong order in one of the likelihood calculations.')

        form.addParam('subtract', BooleanParam, default=True,
                      label='Shift LL matrix by subtracting mean value of each column?',
                      help='This may increase the contrast to help with interpretability. ')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calculationStep')
        self._insertFunctionStep('createOutputStep')


    def calculationStep(self):

        matrix1 = np.load(self.getMatrixPath(1))
        matrix2 = np.load(self.getMatrixPath(2))

        if self.subtract.get():
            matrix1 = np.subtract(matrix1, np.mean(matrix1, axis=0))
            matrix2 = np.subtract(matrix2, np.mean(matrix2, axis=0))

        if self.flipVols.get():
            matrix2 = np.flip(matrix2, axis=0)

        corrcoeffs = np.abs(np.corrcoef(matrix1.flatten(), matrix2.flatten()))[0, 1]
        np.savetxt(self.getCorrCoeffPath(), np.array(corrcoeffs).reshape(-1), fmt='%5.3f')

    def createOutputStep(self):
        # register output files
        self.args = {}

        corrcoeff = np.loadtxt(self.getCorrCoeffPath())

        outSet = EMSet().create(self._getExtraPath())
        item = EMObject()
        setattr(item, '_correlation_coeff', Float(corrcoeff))
        outSet.append(item)

        self.args["outputSet"] = outSet
        self._defineOutputs(**self.args)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        if self.isFinished():
            summary.append("This protocol has finished. The correlation coefficient is %5.3f" \
                           % self.outputSet.getFirstItem()._correlation_coeff)
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            methods.append("This protocol has printed methods")
        return methods

    def _validate(self):
        errors = []
        for i, pointer in enumerate([self.inputParticles1, self.inputParticles2]):
            distanceObject = pointer.get()
            if not isinstance(distanceObject, EMFile):
                if not (hasattr(distanceObject.getFirstItem(), '_xmipp_logLikelihood')
                        or hasattr(distanceObject.getFirstItem(), '_cryolike_logLikelihood')):
                    errors.append('The input particle set {0} must have xmipp or cryolike logLikelihood data'.format(i+1))

        return errors

    def getMatrixPath(self, number):
        if number == 1:
            protocolPath = os.path.dirname(self.inputParticles1.get().getFileName())

        elif number == 2:
            protocolPath = os.path.dirname(self.inputParticles2.get().getFileName())

        return os.path.join(protocolPath, 'extra/matrix.npy')

    def getCorrCoeffPath(self):
        return self._getExtraPath('corrcoeff.txt')
