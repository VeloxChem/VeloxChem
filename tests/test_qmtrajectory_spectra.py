from pathlib import Path

import h5py
import numpy as np
import pytest

from veloxchem.outputstream import OutputStream
from veloxchem.qmtrajectorydriver import QMTrajectoryDriver
from veloxchem.qmtrajectoryanalyzer import QMTrajectoryAnalyzer


@pytest.mark.solvers
def test_qmtrajectorydriver_plot_spectrum_for_lr_and_cpp(monkeypatch):
    plt = pytest.importorskip('matplotlib.pyplot')
    monkeypatch.setattr(plt, 'show', lambda: None)

    drv = QMTrajectoryDriver(ostream=OutputStream(None))

    lr_results = {
        'prop_all': [
            (1, {
                'eigenvalues': np.array([0.20, 0.25]),
                'oscillator_strengths': np.array([0.1, 0.2]),
            }),
            (2, {
                'eigenvalues': np.array([0.21, 0.26]),
                'oscillator_strengths': np.array([0.15, 0.18]),
            }),
        ]
    }
    fig1, ax1 = plt.subplots()
    drv.plot_spectrum(
        lr_results,
        property='absorption',
        x_unit='nm',
        x_range=(300.0, 120.0),
        ax=ax1,
    )
    xlim1 = ax1.get_xlim()
    assert xlim1[0] < xlim1[1]
    avg_line1 = ax1.lines[-1]
    assert np.all(np.diff(avg_line1.get_xdata()) > 0.0)
    legend_labels = [text.get_text() for text in ax1.get_legend().get_texts()]
    assert 'Frames 1-2 (viridis)' in legend_labels
    assert 'Average' in legend_labels
    assert 'States' not in legend_labels
    assert 'Frame 1' not in legend_labels
    assert 'Frame 2' not in legend_labels

    ax1_right = fig1.axes[1]
    y0, _ = ax1_right.get_ylim()
    assert np.isclose(y0, 0.0)

    fig2, ax2 = plt.subplots()
    drv.plot_spectrum(
        lr_results,
        property='absorption',
        x_unit='nm',
        ax=ax2,
    )
    xlim2 = ax2.get_xlim()
    assert xlim2[0] < xlim2[1]
    avg_line2 = ax2.lines[-1]
    assert np.all(np.diff(avg_line2.get_xdata()) > 0.0)

    cpp_results = {
        'prop_all': [
            (1, {
                'frequencies': np.array([0.0, 0.05, 0.10]),
                'response_functions': {
                    ('x', 'x', 0.05): complex(1.0, -0.2),
                    ('y', 'y', 0.05): complex(1.0, -0.2),
                    ('z', 'z', 0.05): complex(1.0, -0.2),
                    ('x', 'x', 0.10): complex(1.0, -0.3),
                    ('y', 'y', 0.10): complex(1.0, -0.3),
                    ('z', 'z', 0.10): complex(1.0, -0.3),
                },
            }),
        ]
    }
    with pytest.raises(AssertionError):
        drv.plot_spectrum(
            cpp_results,
            property='absorption',
            x_unit='nm',
            x_range=(100.0, 250.0),
        )

    fig3, ax3 = plt.subplots()
    drv.plot_spectrum(
        cpp_results,
        property='absorption',
        x_unit='ev',
        x_range=(1.3, 2.8),
        ax=ax3,
    )
    legend_labels_cpp = [text.get_text() for text in ax3.get_legend().get_texts()]
    assert 'Frame 1 (viridis)' in legend_labels_cpp
    assert 'Average' in legend_labels_cpp

    plt.close('all')


@pytest.mark.solvers
def test_qmtrajectoryanalyzer_plot_spectrum_reads_h5(monkeypatch, tmp_path):
    plt = pytest.importorskip('matplotlib.pyplot')
    monkeypatch.setattr(plt, 'show', lambda: None)

    h5_path = tmp_path / 'traj.h5'
    with h5py.File(h5_path, 'w') as hf:
        hf.create_dataset('frame_id', data=np.array([1, 2], dtype=np.int32), maxshape=(None,))
        hf.create_dataset('converged', data=np.array([True, True], dtype=bool), maxshape=(None,))
        hf.create_dataset('basis_set', data=np.array([b'sto-3g', b'sto-3g']))
        hf.create_dataset('dft_func_label', data=np.array([b'b3lyp', b'b3lyp']))
        hf.create_dataset('trajectory_name', data=np.array([b'test']))
        rsp = hf.create_group('rsp')
        rsp.create_dataset('eigenvalues', data=np.array([[0.20, 0.25], [0.21, 0.26]]), maxshape=(None, 2))
        rsp.create_dataset('oscillator_strengths', data=np.array([[0.1, 0.2], [0.15, 0.18]]), maxshape=(None, 2))

    analyzer = QMTrajectoryAnalyzer(h5_path)
    fig, ax = plt.subplots()
    analyzer.plot_spectrum(
        property='absorption',
        x_unit='nm',
        x_range=(300.0, 120.0),
        ax=ax,
    )
    xlim = ax.get_xlim()
    assert xlim[0] < xlim[1]
    avg_line = ax.lines[-1]
    assert np.all(np.diff(avg_line.get_xdata()) > 0.0)

    plt.close('all')
