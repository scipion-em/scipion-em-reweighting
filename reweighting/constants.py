
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

import os
import reweighting

REWEIGHTING_HOME = 'REWEIGHTING_HOME'
REWEIGHTING_URL = 'https://github.com/scipion-em/scipion-em-reweighting'

CONDA_YML = os.path.join(reweighting.__path__[0], 'conda.yaml')

def getReweightingEnvName(version):
    return "reweighting-%s" % version

V0_0_1 = "0.0.1"

VERSIONS = [V0_0_1]
REWEIGHTING_DEFAULT_VER_NUM = V0_0_1

DEFAULT_ENV_NAME = getReweightingEnvName(REWEIGHTING_DEFAULT_VER_NUM)
DEFAULT_ACTIVATION_CMD = 'conda activate ' + DEFAULT_ENV_NAME
REWEIGHTING_ENV_ACTIVATION = 'REWEIGHTING_ENV_ACTIVATION'
