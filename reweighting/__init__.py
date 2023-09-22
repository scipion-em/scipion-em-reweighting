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

import pwem
import reweighting
from reweighting.constants import *
import os
import pyworkflow.utils as pwutils
from scipion.install.funcs import VOID_TGZ

_logo = "icon.png"
_references = ['ensReweighting']

# Use this variable to activate an environment from the Scipion conda
MODEL_CRYOER_ENV_ACTIVATION_VAR = "MODEL_CRYOER_ENV_ACTIVATION"
# Use this general activation variable when installed outside Scipion
MODEL_CRYOER_ACTIVATION_VAR = "MODEL_CRYOER_ACTIVATION"

__version__ = "3.0.1"


class Plugin(pwem.Plugin):
    _homeVar = CRYOER_HOME
    _pathVars = [CRYOER_HOME]
    # We only support our latest release, we can't afford supporting previous releases
    _supportedVersions = [__version__]
    _url = CRYOER_URL

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(MODEL_CRYOER_ACTIVATION_VAR, '')
        cls._defineEmVar(CRYOER_HOME, 'reweighting-0.0.1/CryoER')
        cls._defineVar(MODEL_CRYOER_ENV_ACTIVATION_VAR, cls.getActivationCmd(__version__))

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ)
        return environ

    @classmethod
    def getCRYOERCmd(cls, args):
        cmd = cls.getVar(MODEL_CRYOER_ACTIVATION_VAR)
        if not cmd:
            cmd = cls.getCondaActivationCmd()
            cmd += cls.getVar(MODEL_CRYOER_ENV_ACTIVATION_VAR)
        cmd += " && "
        cmd += args
        return cmd

    @classmethod
    def getActivationCmd(cls, version):
        return 'conda activate {}'.format(cls.getVar(CRYOER_HOME))

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(__version__)

    @classmethod
    def defineBinaries(cls, env):
        def getCondaInstallationCryoER():
            installationCmd = cls.getCondaActivationCmd()
            installationCmd += 'conda env remove -n reweighting && conda env create -n reweighting -f ' + CONDA_YML + " && "
            installationCmd += "conda activate reweighting && "

            clonePath = os.path.join(pwem.Config.EM_ROOT, "torch-batch-svd")
            if not os.path.exists(clonePath):
                installationCmd += "git clone -b master https://github.com/KinglittleQ/torch-batch-svd.git torch-batch-svd && "
            installationCmd += "cd torch-batch-svd && "
            installationCmd += "pip install -Ue . && cd .. && "

            clonePath = os.path.join(pwem.Config.EM_ROOT, "CryoER")
            if not os.path.exists(clonePath):
                installationCmd += "git clone -b jmk_improvements https://github.com/jamesmkrieger/Ensemble-reweighting-using-Cryo-EM-particles.git CryoER && "
            installationCmd += "cd CryoER && "
            installationCmd += "pip install -Ue . && cd .. && "

            installationCmd += "python -c 'import cmdstanpy; cmdstanpy.install_cmdstan()' && "
            
            installationCmd += "touch reweighting_installed"

            return installationCmd

        commands = []
        installationEnv = getCondaInstallationCryoER()
        commands.append((installationEnv, ["reweighting_installed"]))

        env.addPackage('reweighting', version='0.0.1',
                       commands=commands,
                       tar="void.tgz",
                       default=True)
