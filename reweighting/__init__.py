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
import os
import pyworkflow.utils as pwutils
from pyworkflow import Config
from scipion.install.funcs import VOID_TGZ

from reweighting.constants import *

_logo = "icon.png"
_references = ['ensReweighting']
__version__ = "3.0.1"


class Plugin(pwem.Plugin):
    _supportedVersions = VERSIONS
    _url = REWEIGHTING_URL

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(REWEIGHTING_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ)
        return environ

    @classmethod
    def getReweightingCmd(cls, args):
        cmd = '%s %s && ' % (cls.getCondaActivationCmd(), cls.getReweightingEnvActivation())
        cmd += args
        return cmd

    @classmethod
    def getActivationCmd(cls):
        """ Return the activation command. """
        return '%s %s' % (cls.getCondaActivationCmd(),
                          cls.getReweightingEnvActivation())
    
    @classmethod
    def getReweightingEnvActivation(cls):
        """ Activate the conda environment. """
        return cls.getVar(REWEIGHTING_ENV_ACTIVATION)

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(__version__)

    @classmethod
    def defineBinaries(cls, env):
        for ver in VERSIONS:
            cls.addReweightingPackage(env, ver,
                                      default=ver == REWEIGHTING_DEFAULT_VER_NUM)

    @classmethod
    def addReweightingPackage(cls, env, version, default=False):

        def getCondaInstallationReweighting():
            ENV_NAME = getReweightingEnvName(version)
            installationCmd = cls.getCondaActivationCmd()
            installationCmd += f" conda env create -n {ENV_NAME} -f {CONDA_YML} --force && "
            installationCmd += f"conda activate {ENV_NAME} && "

            clonePath = os.path.join(pwem.Config.EM_ROOT, "Reweighting")
            if not os.path.exists(clonePath):
                installationCmd += "git clone -b main https://github.com/jamesmkrieger/Ensemble-reweighting-using-Cryo-EM-particles.git Reweighting && "

            installationCmd += "cd Reweighting && "
            installationCmd += "pip install -Ue . && cd .. && "
            installationCmd += "python -c 'import cmdstanpy; cmdstanpy.install_cmdstan()' && "

            installationCmd += "touch reweighting_installed"
            return installationCmd

        def getCondaInstallationTorchSVD():
            ENV_NAME = getReweightingEnvName(version)
            installationCmd = cls.getCondaActivationCmd()
            installationCmd += f"conda activate {ENV_NAME} && "

            clonePath = os.path.join(pwem.Config.EM_ROOT, "torch-batch-svd")
            if not os.path.exists(clonePath):
                installationCmd += "git clone -b master https://github.com/KinglittleQ/torch-batch-svd.git torch-batch-svd && "
            installationCmd += "cd torch-batch-svd && "
            installationCmd += "pip install -Ue . && cd .. && "

            installationCmd += "touch reweighting_torch_svd_installed"
            return installationCmd

        commands = []
        installationEnv = getCondaInstallationReweighting()
        installationTorchSVD = getCondaInstallationTorchSVD()
        commands.append((installationEnv, ["reweighting_installed"]))
        commands.append((installationTorchSVD, ["reweighting_torch_svd_installed"]))

        env.addPackage('reweighting', version=version,
                       commands=commands,
                       tar="void.tgz",
                       default=True)
