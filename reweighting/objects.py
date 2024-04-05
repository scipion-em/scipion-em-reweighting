# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
from pwem.objects import CTFModel, Image
from pyworkflow.object import String

class ReweightingCTF(CTFModel):

    def __init__(self, **kwargs):
        CTFModel.__init__(self, **kwargs)
        self._ctfFile = String()
   
    def copyInfo(self, other):
        self.copyAttributes(other, '_defocusU', '_defocusV', '_defocusAngle',
                            '_defocusRatio', '_psdFile', '_micFile',
                            '_resolution', '_fitQuality', '_ctfFile')
        if other.hasPhaseShift():
            self.setPhaseShift(other.getPhaseShift())   

    def getCtfFile(self):
        return self._ctfFile.get()

    def setCtfFile(self, value):
        self._ctfFile.set(value)

    def calcCtfImage(self):
        """ Calculate CTF image from existing data in object.

        This could be Fourier transform of _psdFile data or 
        fresh calculation from defocus, amp contrast, b_factor, etc.
        """
        pass