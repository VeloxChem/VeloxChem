from pathlib import Path

import pytest

from veloxchem.outputstream import OutputStream
from veloxchem.qmtrajectoryparser import QMTrajectoryParser


@pytest.mark.solvers
def test_qmtrajectoryparser_stride_on_xyz(tmp_path):
    xyz_path = tmp_path / 'traj.xyz'
    xyz_path.write_text(
        '1\nframe 0\nH 0.0 0.0 0.0\n'
        '1\nframe 1\nH 0.0 0.0 1.0\n'
        '1\nframe 2\nH 0.0 0.0 2.0\n'
        '1\nframe 3\nH 0.0 0.0 3.0\n'
        '1\nframe 4\nH 0.0 0.0 4.0\n'
        '1\nframe 5\nH 0.0 0.0 5.0\n'
    )

    parser = QMTrajectoryParser(ostream=OutputStream(None))
    frames = parser.structures(
        trajectory_file=str(xyz_path),
        qm_charge=0,
        qm_multiplicity=1,
        start_frame=1,
        end_frame=6,
        stride=2,
    )

    assert [frame['frame'] for frame in frames] == [1, 3, 5]
    assert [frame['frame_index'] for frame in frames] == [0, 1, 2]
    assert all(frame['total_frames'] == 3 for frame in frames)


@pytest.mark.solvers
def test_qmtrajectoryparser_invalid_stride_on_xyz(tmp_path):
    xyz_path = tmp_path / 'traj.xyz'
    xyz_path.write_text(
        '1\nframe 0\nH 0.0 0.0 0.0\n'
        '1\nframe 1\nH 0.0 0.0 1.0\n'
    )

    parser = QMTrajectoryParser(ostream=OutputStream(None))

    with pytest.raises(ValueError, match='stride'):
        parser.structures(
            trajectory_file=str(xyz_path),
            qm_charge=0,
            qm_multiplicity=1,
            stride=0,
        )
