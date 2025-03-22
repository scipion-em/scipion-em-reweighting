# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
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
This module implements visualization for scatter plots to correlate likelihoods.
"""
import matplotlib.pyplot as plt
import numpy as np

from pwem.viewers.plotter import EmPlotter

from pyworkflow.protocol.params import LabelParam, IntParam, FloatParam, BooleanParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from reweighting.protocols.protocol_correlate_likelihoods import ReweightingCorrelateProtocol

_invalidInputStr = 'Invalid input'

class ReweightingCorrelationViewer(ProtocolViewer):
    """ Visualization of results from the Xmipp log likelihood protocol.
    """
    _label = 'Correlation viewer'
    _targets = [ReweightingCorrelateProtocol]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')

        group = form.addGroup('Volumes range')
        group.addParam('volNumber1', IntParam, default=-1,
                      label='Initial volume number')
        group.addParam('volNumber2', IntParam, default=-1,
                      label='Final volume number')

        group.addParam('label', BooleanParam, label="Show correlation coefficient?", default=True,
                      help='Select whether to show the Pearson correlation coefficient.')

        group.addParam('flipVols', BooleanParam, label="Flip volumes for matrix 2?", default=False,
                      help='Select whether to flip volumes for matrix 2. This may be useful '
                           'if they are in the wrong order in one of the likelihood calculations.')

        group.addParam('subtract', BooleanParam, default=True,
                      label='Shift LL matrix by subtracting mean value of each column?',
                      help='This may increase the contrast to help with interpretability. ')

        form.addParam('displayScatterPlot', LabelParam, default=False,
                label="Plot scatter plot to illustrate correlations?",
                help="Scatter plots use flattened matrices from the selected volumes.")
        
    def _getVisualizeDict(self):
        return {'displayScatterPlot': self._viewScatter} 

    def _viewScatter(self, paramName):
        """ visualization of scatters and correlation coefficients for all references or the range selected. """


        self._checkNumbers(1)
        matrix1 = np.load(self.protocol.getMatrixPath(1))[self.volumeNumber1:self.volumeNumber2]

        self._checkNumbers(2)
        matrix2 = np.load(self.protocol.getMatrixPath(2))[self.volumeNumber1:self.volumeNumber2]

        if self.flipVols.get():
            matrix2 = np.flip(matrix2, axis=0)

        if self.subtract.get():
            matrix1 = np.subtract(matrix1, np.mean(matrix1, axis=0))
            matrix2 = np.subtract(matrix2, np.mean(matrix2, axis=0))

        plotter = EmPlotter()
        corrcoeff = np.corrcoef(matrix1.flatten(), matrix2.flatten())[0,1]
        plt.scatter(matrix1.flatten(), matrix2.flatten(), label='%6.3f' % corrcoeff)
        if self.label.get():
            plt.legend()

        plt.xlabel('Likelihood 1')
        plt.ylabel('Likelihood 2')

        return [plotter]


    def _checkNumbers(self, setNumber):

        self.volumeNumber1 = self.volNumber1.get()-1 if self.volNumber1.get() != -1 else 0

        if setNumber == 1:
            string = 'particle set 1'
            items = self.protocol.inputParticles1.get()
        elif setNumber == 2:
            string = 'particle set 2'
            items = self.protocol.inputParticles2.get()

        self.volumeNumber2 = self.volNumber2.get() if self.volNumber1.get() != -1 else len(items)+1

        if self.volumeNumber1+1 > self.volumeNumber1:
            return [self.errorMessage("Invalid {0} range\n"
                                      "Initial {0} number can not be " 
                                      "bigger than the final one.".format(string), 
                                      title=_invalidInputStr)]

        elif self.volumeNumber1 < -1:
            return [self.errorMessage("Invalid {0} range\n"
                                      "Initial {0} number can not be " 
                                      "smaller than -1.".format(string), 
                                      title=_invalidInputStr)]

        elif self.volumeNumber1 < -1:
            return [self.errorMessage("Invalid {0} range\n"
                                      "Final {0} number can not be " 
                                      "smaller than -1.".format(string), 
                                      title=_invalidInputStr)]
        
        if self.volumeNumber1 != -1:
            try:
                _ = items[self.volumeNumber1]
            except IndexError:
                return [self.errorMessage("Invalid initial {0} number {1}\n"
                                         "Display the output {0}s to see "
                                         "the availables ones.".format(string, self.volNumber1.get()+1),
                                         title=_invalidInputStr)]
            
        if self.volumeNumber2 != -1:
            try:
                _ = items[self.volumeNumber2-1]
            except IndexError:
                return [self.errorMessage("Invalid final {0} number {1}\n"
                                         "Display the output {0}s to see "
                                         "the availables ones.".format(string, self.volNumber2.get()),
                                         title=_invalidInputStr)]
