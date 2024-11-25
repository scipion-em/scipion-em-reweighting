# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Slavica Jonic  (slavica.jonic@upmc.fr)
# *              James Krieger (jmkrieger@cnb.csic.es)
# *              Ricardo Serrano Gutiérrez (rserranogut@hotmail.com)  
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
This module implement the wrappers around ProDy GNM 
visualization programs.
"""
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os

from pwem.viewers.plotter import EmPlotter

from pyworkflow.protocol.params import LabelParam, IntParam, FloatParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from reweighting.protocols.protocol_estimate import ReweightingEstimateWeightsProtocol

_invalidInputStr = 'Invalid input'

class ReweightingWeightsViewer(ProtocolViewer):
    """ Visualization of results from the Xmipp log likelihood protocol.
    """
    _label = 'Log likelihood matrix viewer'
    _targets = [ReweightingEstimateWeightsProtocol]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def _defineParams(self, form):
        self.outputs = self.protocol.outputSet

        form.addSection(label='Visualization')

        group = form.addGroup('Volumes range')
        group.addParam('volNumber1', IntParam, default=-1,
                      label='Initial volume number')
        group.addParam('volNumber2', IntParam, default=-1,
                      label='Final volume number')

        group = form.addGroup('Values range')
        group.addParam('vmin', FloatParam, default=-1,
                      label='Minimum value',
                      help='The axis will be cropped below this value')
        group.addParam('vmax', FloatParam, default=-1,
                      label='Maximum value',
                      help='The axis will be cropped above this value')

        form.addParam('displayWeights', LabelParam, default=False,
                label="Plot weight histogram?",
                help="Bar plots with errors represent mean and std weights.")
        
    def _getVisualizeDict(self):
        return {'displayWeights': self._viewWeights} 

    def _viewWeights(self, paramName):
        """ visualization of mean and std weights for all references or the range selected. """
        volNumber1 = self.volNumber1.get()-1 if self.volNumber1.get() != -1 else -1
        volNumber2 = self.volNumber2.get() # no subtraction as end of range
        self._checkNumbers(volNumber1, volNumber2, 'volume')
        
        means = [item._reweightingMean.get() for item in self.outputs]
        stds = [item._reweightingStd.get() for item in self.outputs]

        if volNumber1 != -1 and volNumber2 != -1:
            means = means[volNumber1:volNumber2]
            stds = stds[volNumber1:volNumber2]
            volMin = volNumber1
            volMax = volNumber2
        elif volNumber1 != -1:
            means = means[volNumber1:]
            stds = stds[volNumber1:]
            volMin = volNumber1
            volMax = len(self.outputs)+1
        elif volNumber2 != -1:
            means = means[:volNumber2]
            stds = stds[:volNumber2]
            volMin = 1
            volMax = volNumber2
        else:
            volMin = 1
            volMax = len(self.outputs)+1

        vmin = self.vmin.get()
        if vmin == -1:
            vmin = None

        vmax = self.vmax.get()
        if vmax == -1:
            vmax = None

        plotter = EmPlotter()
        plt.bar(range(len(means)), means)
        plt.errorbar(range(len(means)), means, yerr=stds, 
                     fmt='none', color='k', capsize=5)
        
        labels = ['{:5.3f} +/- {:6.4f}'.format(m, s) for (m, s) in zip(means, stds)]
        for i in range(len(means)):
            plt.text(i, means[i]+0.05, labels[i],
                     horizontalalignment='center')

        plt.xlabel('Reference volumes')
        plt.xticks(range(volMax-volMin), range(volMin, volMax))

        plt.ylabel('Weights')
        plt.ylim([vmin, vmax])

        return [plotter]


    def _checkNumbers(self, number1, number2, string):

        items = self.outputs

        if number1+1 > number2:
            return [self.errorMessage("Invalid {0} range\n"
                                      "Initial {0} number can not be " 
                                      "bigger than the final one.".format(string), 
                                      title=_invalidInputStr)]

        elif number1 < -1:
            return [self.errorMessage("Invalid {0} range\n"
                                      "Initial {0} number can not be " 
                                      "smaller than -1.".format(string), 
                                      title=_invalidInputStr)]

        elif number2 < -1:
            return [self.errorMessage("Invalid {0} range\n"
                                      "Final {0} number can not be " 
                                      "smaller than -1.".format(string), 
                                      title=_invalidInputStr)]
        
        if number1 != -1:
            try:
                _ = items[number1]
            except IndexError:
                return [self.errorMessage("Invalid initial {0} number {1}\n"
                                         "Display the output {0}s to see "
                                         "the availables ones.".format(string, number1+1),
                                         title=_invalidInputStr)]
            
        if number2 != -1:
            try:
                _ = items[number2-1]
            except IndexError:
                return [self.errorMessage("Invalid final {0} number {1}\n"
                                         "Display the output {0}s to see "
                                         "the availables ones.".format(string, number2),
                                         title=_invalidInputStr)]
