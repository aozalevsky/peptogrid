# This file is part of the PeptoGrid package for the rescoring
# of AutoDock Vina docking poses
#
# Copyright (c) 2017-2018, by Arthur Zalevsky <aozalevsky@fbb.msu.ru>
#
# AffBio is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# AffBio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with AffBio; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#

import argparse as ag
import ConfigParser
import numpy as np


np.set_printoptions(threshold='nan')


def get_bounding(R, padding, step):
    minXYZ = np.amin(R, axis=0)
    maxXYZ = np.amax(R, axis=0)
    minXYZ -= (padding + step)
    maxXYZ += (padding + step)
    shape = np.ceil((maxXYZ - minXYZ) / step).astype(np.int)
    minXYZ = adjust_grid(minXYZ, step)
    return(minXYZ, shape)


def adjust_grid(c, step, padding=0):
    nc = c - padding
    nc /= step
    nc = np.floor(nc) * step
    return nc


def sphere2grid(c, r, step, value=-np.inf):
    nc = adjust_grid(c, step)
    adj = nc - c
    r2x = np.arange(-r - adj[0], r - adj[0] + step, step) ** 2
    r2y = np.arange(-r - adj[1], r - adj[1] + step, step) ** 2
    r2z = np.arange(-r - adj[2], r - adj[2] + step, step) ** 2
    dist2 = r2x[:, None, None] + r2y[:, None] + r2z
    nc = nc - dist2.shape[0] / 2 * step
    vol = (dist2 <= r ** 2).astype(np.int)
    vol *= value
    return (nc, vol)


#  http://www.weslack.com/question/1552900000001895550

def process_atom_oddt(A, step):

    AminXYZ, Agrid = sphere2grid(
        A[1], avdw[A[5][0]], step, 1)

    np.clip(Agrid, 0, 1, out=Agrid)
    return (Agrid, AminXYZ)

avdw = {
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'P': 1.8,
    'S': 1.8,
}

vina_types = (
    'C.2',
    'C.3',
    'C.ar',
    'C.cat',
    'O.2',
    'O.3',
    'O.co2',
    'N.3',
    'N.4',
    'N.ar',
    'N.am',
    'N.pl3',
    'S.3',
)


# http://stackoverflow.com/questions/2819696/parsing-properties-file-in-python/2819788#2819788

class FakeSecHead(object):

    def __init__(self, fp):
        self.fp = fp
        self.sechead = '[asection]\n'

    def readline(self):
        if self.sechead:
            try:
                return self.sechead
            finally:
                self.sechead = None
        else:
            return self.fp.readline()


def get_box_from_vina(f):
    cp = ConfigParser.SafeConfigParser()
    cp.readfp(FakeSecHead(open(f)))
    cfg = dict(cp.items('asection'))

    center = np.array([
                      cfg['center_x'],
                      cfg['center_y'],
                      cfg['center_z'],
                      ],
                      dtype=np.float)

    box = np.array([
                   cfg['size_x'],
                   cfg['size_y'],
                   cfg['size_z'],
                   ],
                   dtype=np.float)

    return (center, box)

one_let = [
    "A", "R", "N", "D", "C", "E", "Q", "G", "H",
    "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "BCK", "UNK"]


three_let = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "BCK", "UNK", "NME", "ACE",
]


def get_args_grid():
    """Parse cli arguments"""

    parser = ag.ArgumentParser(
        description='Grid scripts')

    parser.add_argument('-m',
                        required=True,
                        dest='Sfn',
                        metavar='OUTPUT.hdf5',
                        help='Output grid file stored in HDF5 matrix')

    parser.add_argument('-f',
                        nargs='+',
                        type=str,
                        dest='pdb_list',
                        metavar='FILE',
                        help='PDBQT files')

    parser.add_argument('-s', '--step',
                        type=float,
                        dest='step',
                        metavar='STEP',
                        default=0.1,
                        help='Grid step in Angstroms')

    parser.add_argument('-p', '--padding',
                        type=int,
                        dest='pad',
                        metavar='PADDING',
                        default=1,
                        help='Padding around cell in Angstroms')

    parser.add_argument('-c', '--config',
                        type=str,
                        dest='config',
                        metavar='CONFIG',
                        help='Vina config file')

    args = parser.parse_args()

    args_dict = vars(args)

    return args_dict


def get_args_score():
    """Parse cli arguments"""
    parser = ag.ArgumentParser(
        description='Grid scripts')

    parser.add_argument('-m',
                        required=True,
                        dest='Sfn',
                        metavar='FILE.hdf5',
                        help='PeptoGrid HDF5 file for all matrices')

    parser.add_argument('-o', '--output',
                        dest='output',
                        metavar='OUTPUT',
                        help='Output score file name',
                        )

    parser.add_argument('-f',
                        nargs='+',
                        type=str,
                        dest='pdb_list',
                        metavar='FILE',
                        help='PDBQT files to score')

    parser.add_argument('-e', '--exclude',
                        nargs='*',
                        type=str,
                        dest='excl',
                        metavar='TYPE',
                        help='Vina types to exclude from scoring')

    parser.add_argument('-i', '--include',
                        nargs='*',
                        type=str,
                        dest='incl',
                        metavar='TYPE',
                        help='Vina types to exclude from scoring')

    parser.add_argument('-s', '--step',
                        type=float,
                        dest='step',
                        metavar='STEP',
                        default=0.1,
                        help='Grid step in Angstroms')

    parser.add_argument('-p', '--padding',
                        type=int,
                        dest='pad',
                        metavar='PADDING',
                        default=1,
                        help='Padding around cell in Angstroms')

    parser.add_argument('-c', '--config',
                        type=str,
                        dest='config',
                        metavar='CONFIG',
                        help='Vina config file')

    args = parser.parse_args()

    args_dict = vars(args)

    return args_dict
