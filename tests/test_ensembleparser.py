import pytest
import numpy as np

import veloxchem.ensembleparser as ensembleparser_module
from veloxchem.ensembleparser import EnsembleParser

pytest.importorskip("MDAnalysis")


def make_simple_tip3_pdb(tmp_path):
    """Writes a simple PDB containing one QM atom and one TIP3 water."""

    pdb_path = tmp_path / "simple_tip3.pdb"
    pdb_lines = [
        _make_pdb_record(1, "C1", "LIG", "A", 1, 0.0, 0.0, 0.0, "C", record_name="HETATM"),
        _make_pdb_record(2, "OH2", "TIP3", "A", 2, 2.0, 0.0, 0.0, "O", record_name="HETATM"),
        _make_pdb_record(3, "H1", "TIP3", "A", 2, 2.7, 0.0, 0.0, "H", record_name="HETATM"),
        _make_pdb_record(4, "H2", "TIP3", "A", 2, 1.7, 0.7, 0.0, "H", record_name="HETATM"),
        "END",
    ]
    pdb_path.write_text("\n".join(pdb_lines) + "\n")
    return pdb_path


def _make_pdb_record(serial, atom_name, residue_name, chain, resid, x, y, z, element,
                     record_name="ATOM"):
    """Format a single ATOM/HETATM record for a PDB file."""
    return (
        f"{record_name:<6s}{serial:5d} "
        f"{atom_name:^4s} "
        f"{residue_name:>4s}"
        f"{chain:1s}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}          "
        f"{element:>2s}"
    )


class TestEnsembleParser:

    @pytest.mark.parametrize(
        ("mda_version", "expected_kwargs"),
        [
            ("2.7.9", {"guess_bonds": True}),
            (
                "2.8.0",
                {"to_guess": ("bonds", "angles", "dihedrals")},
            ),
        ],
    )
    def test_mdanalysis_pdb_version_options(
        self, monkeypatch, mda_version, expected_kwargs
    ):
        class UniverseCalled(Exception):
            pass

        universe_kwargs = {}

        def fake_version(distribution_name):
            assert distribution_name == "MDAnalysis"
            return mda_version

        def fake_universe(*args, **kwargs):
            universe_kwargs.update(kwargs)
            raise UniverseCalled

        monkeypatch.setattr(ensembleparser_module, "version", fake_version)
        monkeypatch.setattr(
            ensembleparser_module.mda,
            "Universe",
            fake_universe,
        )

        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        with pytest.raises(UniverseCalled):
            ens_parser.structures(
                pdb_file="structure.pdb",
                qm_region="resname LIG",
            )

        assert universe_kwargs == expected_kwargs

    def test_pdb_file_and_tip3_pe_pot_files(self, tmp_path):
        from veloxchem.veloxchemlib import mpi_master
        from veloxchem.ensembledriver import EnsembleDriver

        simple_tip3_pdb = make_simple_tip3_pdb(tmp_path)
        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        snapshots = ens_parser.structures(
            pdb_file=simple_tip3_pdb,
            qm_region="resname LIG",
            pe_cutoff=3.0,
        )

        assert len(snapshots) == 1
        assert snapshots[0]["qm_atom_names"].tolist() == ["C1"]
        assert snapshots[0]["pe_atom_names"].tolist() == ["OH2", "H1", "H2"]
        assert snapshots[0]["pe_resnames"].tolist() == ["TIP3"] * 3

        tip3p_snapshot = dict(snapshots[0])
        tip3p_snapshot["frame"] = 1
        tip3p_snapshot["pe_resnames"] = np.array(["TIP3P"] * 3, dtype=object)

        ens_drv = EnsembleDriver()
        ens_drv.ostream.mute()
        ens_drv.set_env_models(pe_model="SEP")
        ens_drv.write_pot_files(
            [snapshots[0], tip3p_snapshot],
            outdir=tmp_path,
        )

        if ens_drv.rank == mpi_master():
            for frame, resname in ((0, "TIP3"), (1, "TIP3P")):
                pot_text = (
                    tmp_path / f"pe_frame_{frame:06d}.pot"
                ).read_text()
                charge_lines = pot_text.split(
                    "@charges\n", 1
                )[1].split("@end", 1)[0].splitlines()

                assert [line.split() for line in charge_lines] == [
                    ["O", "-0.67444000", f"{resname}_pe"],
                    ["H", "0.33722000", f"{resname}_pe"],
                    ["H", "0.33722000", f"{resname}_pe"],
                ]
                assert "@polarizabilities\n" in pot_text

    def test_tip3_npe_only_pot_files(self, tmp_path):
        from veloxchem.veloxchemlib import mpi_master
        from veloxchem.ensembledriver import EnsembleDriver

        simple_tip3_pdb = make_simple_tip3_pdb(tmp_path)
        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        snapshots = ens_parser.structures(
            pdb_file=simple_tip3_pdb,
            qm_region="resname LIG",
            npe_cutoff=3.0,
        )

        tip3p_snapshot = dict(snapshots[0])
        tip3p_snapshot["frame"] = 1
        tip3p_snapshot["npe_resnames"] = np.array(["TIP3P"] * 3, dtype=object)

        ens_drv = EnsembleDriver()
        ens_drv.ostream.mute()
        ens_drv.set_env_models(pe_model="SEP", npe_model="tip3p")
        ens_drv.write_pot_files(
            [snapshots[0], tip3p_snapshot],
            outdir=tmp_path,
        )

        if ens_drv.rank == mpi_master():
            for frame, resname in ((0, "TIP3"), (1, "TIP3P")):
                pot_text = (
                    tmp_path / f"pe_frame_{frame:06d}.pot"
                ).read_text()
                charge_lines = pot_text.split(
                    "@charges\n", 1
                )[1].split("@end", 1)[0].splitlines()

                assert [line.split() for line in charge_lines] == [
                    ["O", "-0.83400000", f"{resname}_npe"],
                    ["H", "0.41700000", f"{resname}_npe"],
                    ["H", "0.41700000", f"{resname}_npe"],
                ]
                assert "@polarizabilities\n" not in pot_text

    # _histidine_resname_from_atoms (static method)

    @pytest.mark.parametrize(
        ("resname", "atom_names", "expected"),
        [
            # delta-protonated (HD1 only)
            ("HIS", {"HD1"}, "HSD"),
            ("HIS", {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "HD1"}, "HSD"),
            # epsilon-protonated (HE2 only)
            ("HIS", {"HE2"}, "HSE"),
            ("HIS", {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "HE2"}, "HSE"),
            # epsilon-protonated with alternative HNE2 name
            ("HIS", {"HNE2"}, "HSE"),
            ("HIS", {"HD1", "HNE2"}, "HSP"),
            # doubly protonated (both HD1 and HE2)
            ("HIS", {"HD1", "HE2"}, "HSP"),
            ("HIS", {"N", "CA", "C", "O", "HD1", "HE2"}, "HSP"),
            # ambiguous (neither HD1 nor HE2)
            ("HIS", {"N", "CA", "C", "O"}, None),
            ("HIS", set(), None),
            # N-terminal prefixed
            ("NHIS", {"HD1"}, "NHSD"),
            ("NHIS", {"HE2"}, "NHSE"),
            ("NHIS", {"HD1", "HE2"}, "NHSP"),
            ("NHIS", {"N", "CA", "C"}, None),
            # C-terminal prefixed
            ("CHIS", {"HD1"}, "CHSD"),
            ("CHIS", {"HE2"}, "CHSE"),
            ("CHIS", {"HD1", "HE2"}, "CHSP"),
            # non-HIS residues (should not be renamed)
            ("GLU", {"HD1"}, None),
            ("ASP", {"HE2"}, None),
            ("CYS", {"HD1", "HE2"}, None),
            ("ALA", {"HD1"}, None),
            # already-specific CHARMM names (core != "HIS")
            ("HSD", {"HD1"}, None),
            ("HSE", {"HE2"}, None),
            ("HSP", {"HD1", "HE2"}, None),
            # defensive: non-string resname
            ("HIS", frozenset(["HD1"]), "HSD"),
        ],
    )
    def test_histidine_resname_from_atoms(self, resname, atom_names, expected):
        """Static method infers HSD/HSE/HSP from atom names."""
        result = EnsembleParser._histidine_resname_from_atoms(resname, atom_names)
        assert result == expected

    # _protonation_resname_map with a multi-residue PDB

    @staticmethod
    def _make_protonation_pdb(tmp_path):
        """Write a PDB covering protonation-state residue renaming."""

        records = []
        serial = 0
        residue_atom_counts = {}

        def atom(name, resname, resid, element="C"):
            nonlocal serial
            serial += 1
            atom_index = residue_atom_counts.get(resid, 0)
            residue_atom_counts[resid] = atom_index + 1
            # Use residue-separated offsets to keep generated atoms distinct.
            x = float(resid) * 20.0 + float(atom_index) * 1.2
            records.append(_make_pdb_record(serial, name, resname, "A", resid, x, 0.0, 0.0, element))

        # Residue 1: HIS with HD1 -> HSD
        atom("N", "HIS", 1, "N")
        atom("CA", "HIS", 1, "C")
        atom("HD1", "HIS", 1, "H")

        # Residue 2: HIS with HE2 -> HSE
        atom("N", "HIS", 2, "N")
        atom("CA", "HIS", 2, "C")
        atom("HE2", "HIS", 2, "H")

        # Residue 3: HIS with HNE2 (alternative name) -> HSE
        atom("N", "HIS", 3, "N")
        atom("CA", "HIS", 3, "C")
        atom("HNE2", "HIS", 3, "H")

        # Residue 4: HIS with both HD1 and HE2 -> HSP
        atom("N", "HIS", 4, "N")
        atom("HD1", "HIS", 4, "H")
        atom("HE2", "HIS", 4, "H")

        # Residue 5: HIS without HD1/HE2 -> unmapped
        atom("N", "HIS", 5, "N")
        atom("CA", "HIS", 5, "C")

        # Residue 6: GLU with HE1/HE2 -> GLH
        atom("N", "GLU", 6, "N")
        atom("HE1", "GLU", 6, "H")
        atom("HE2", "GLU", 6, "H")

        # Residue 7: ASP with HD2 -> ASH
        atom("N", "ASP", 7, "N")
        atom("HD2", "ASP", 7, "H")

        # Residue 8: CYS without HG/HG1 -> CYX
        atom("N", "CYS", 8, "N")
        atom("SG", "CYS", 8, "S")

        # Residue 9: CYS with HG -> normal CYS (unmapped)
        atom("N", "CYS", 9, "N")
        atom("HG", "CYS", 9, "H")

        records.append("END")
        pdb_path = tmp_path / "protonation.pdb"
        pdb_path.write_text("\n".join(records) + "\n")
        return pdb_path

    def test_protonation_resname_map(self, tmp_path):
        """Full _protonation_resname_map integration covering HIS/GLU/ASP/CYS."""
        import MDAnalysis as mda

        pdb_path = self._make_protonation_pdb(tmp_path)
        u = mda.Universe(str(pdb_path))
        env_atoms = u.select_atoms("all")

        ens_parser = EnsembleParser()
        prot_map = ens_parser._protonation_resname_map(env_atoms)

        assert prot_map[0] == "HSD"
        assert prot_map[1] == "HSE"
        assert prot_map[2] == "HSE"   # HNE2 → HSE
        assert prot_map[3] == "HSP"
        # Residue 4 (resindex=4): HIS without HD1/HE2 → not in map
        assert 4 not in prot_map
        # Existing rules still work
        assert prot_map[5] == "GLH"
        assert prot_map[6] == "ASH"
        assert prot_map[7] == "CYX"
        # CYS with HG stays unmapped
        assert 8 not in prot_map

        # No unexpected entries
        assert len(prot_map) == 7

    # _prefixed_resname composition with protonation-corrected names

    def test_prefixed_resname_with_protonation_names(self):
        """_prefixed_resname composes terminal prefixes with protonation names."""
        ens_parser = EnsembleParser()

        # N-terminal prefix + protonation-corrected names
        assert ens_parser._prefixed_resname("HSD", "N") == "NHSD"
        assert ens_parser._prefixed_resname("HSE", "N") == "NHSE"
        assert ens_parser._prefixed_resname("HSP", "N") == "NHSP"

        # C-terminal prefix + protonation-corrected names
        assert ens_parser._prefixed_resname("HSD", "C") == "CHSD"
        assert ens_parser._prefixed_resname("HSE", "C") == "CHSE"
        assert ens_parser._prefixed_resname("HSP", "C") == "CHSP"

        # Already-prefixed names are returned unchanged (idempotent)
        assert ens_parser._prefixed_resname("NHSD", "N") == "NHSD"
        assert ens_parser._prefixed_resname("CHSE", "C") == "CHSE"
        assert ens_parser._prefixed_resname("NHSP", "C") == "NHSP"

        # Other protonatable residues also compose correctly
        assert ens_parser._prefixed_resname("GLH", "N") == "NGLH"
        assert ens_parser._prefixed_resname("ASH", "C") == "CASH"
        assert ens_parser._prefixed_resname("CYX", "N") == "NCYX"

    # Terminal + protonation composition via PDB

    @staticmethod
    def _make_terminal_protonation_pdb(tmp_path):
        """Write a PDB covering combined terminal and protonation renaming."""

        records = []
        serial = 0
        chain_atom_counts = {}

        def atom(name, resname, chain, resid, element="C"):
            nonlocal serial
            serial += 1
            atom_index = chain_atom_counts.get(chain, 0)
            chain_atom_counts[chain] = atom_index + 1
            # Use floating-point offsets to keep generated atoms distinct by chain.
            x = float(ord(chain) - ord("A")) * 20.0 + float(atom_index) * 1.2
            records.append(_make_pdb_record(serial, name, resname, chain, resid, x, 0.0, 0.0, element))

        # Chain A: N-terminal HIS (H1,H2,H3) + HD1  and  ALA spacer
        atom("N", "HIS", "A", 1, "N")
        atom("H1", "HIS", "A", 1, "H")
        atom("H2", "HIS", "A", 1, "H")
        atom("H3", "HIS", "A", 1, "H")
        atom("CA", "HIS", "A", 1, "C")
        atom("HD1", "HIS", "A", 1, "H")
        atom("N", "ALA", "A", 2, "N")
        atom("CA", "ALA", "A", 2, "C")

        # Chain B: C-terminal HIS (OT1) + HE2  and  ALA spacer
        atom("N", "ALA", "B", 1, "N")
        atom("CA", "ALA", "B", 1, "C")
        atom("N", "HIS", "B", 2, "N")
        atom("CA", "HIS", "B", 2, "C")
        atom("C", "HIS", "B", 2, "C")
        atom("O", "HIS", "B", 2, "O")
        atom("OT1", "HIS", "B", 2, "O")
        atom("HE2", "HIS", "B", 2, "H")

        # Chain C: N-terminal HIS with both HD1+HE2  and  ALA spacer
        atom("N", "HIS", "C", 1, "N")
        atom("H1", "HIS", "C", 1, "H")
        atom("H2", "HIS", "C", 1, "H")
        atom("H3", "HIS", "C", 1, "H")
        atom("CA", "HIS", "C", 1, "C")
        atom("HD1", "HIS", "C", 1, "H")
        atom("HE2", "HIS", "C", 1, "H")
        atom("N", "ALA", "C", 2, "N")
        atom("CA", "ALA", "C", 2, "C")

        # Chain D: N-terminal HIS without HD1/HE2  and  ALA spacer
        atom("N", "HIS", "D", 1, "N")
        atom("H1", "HIS", "D", 1, "H")
        atom("H2", "HIS", "D", 1, "H")
        atom("H3", "HIS", "D", 1, "H")
        atom("CA", "HIS", "D", 1, "C")
        atom("N", "ALA", "D", 2, "N")
        atom("CA", "ALA", "D", 2, "C")

        # Chain E: N-terminal GLU with HE1/HE2  and  ALA spacer
        atom("N", "GLU", "E", 1, "N")
        atom("H1", "GLU", "E", 1, "H")
        atom("H2", "GLU", "E", 1, "H")
        atom("H3", "GLU", "E", 1, "H")
        atom("CA", "GLU", "E", 1, "C")
        atom("HE1", "GLU", "E", 1, "H")
        atom("HE2", "GLU", "E", 1, "H")
        atom("N", "ALA", "E", 2, "N")
        atom("CA", "ALA", "E", 2, "C")

        # Chain F: C-terminal ASP (OT1) + HD2  and  ALA spacer
        atom("N", "ALA", "F", 1, "N")
        atom("CA", "ALA", "F", 1, "C")
        atom("N", "ASP", "F", 2, "N")
        atom("CA", "ASP", "F", 2, "C")
        atom("C", "ASP", "F", 2, "C")
        atom("O", "ASP", "F", 2, "O")
        atom("OT1", "ASP", "F", 2, "O")
        atom("HD2", "ASP", "F", 2, "H")

        # Chain G: N-terminal CYS without HG/HG1  and  ALA spacer
        atom("N", "CYS", "G", 1, "N")
        atom("H1", "CYS", "G", 1, "H")
        atom("H2", "CYS", "G", 1, "H")
        atom("H3", "CYS", "G", 1, "H")
        atom("CA", "CYS", "G", 1, "C")
        atom("SG", "CYS", "G", 1, "S")
        atom("N", "ALA", "G", 2, "N")
        atom("CA", "ALA", "G", 2, "C")

        for name, x, y, z, element in (
            ("OH2", 0.0, 5.0, 0.0, "O"),
            ("H1", 0.7, 5.0, 0.0, "H"),
            ("H2", -0.3, 5.7, 0.0, "H"),
        ):
            serial += 1
            records.append(_make_pdb_record(
                serial, name, "TIP3", "Z", 999, x, y, z, element,
                record_name="HETATM"))

        records.append("END")
        pdb_path = tmp_path / "terminal_protonation.pdb"
        pdb_path.write_text("\n".join(records) + "\n")
        return pdb_path

    def test_terminal_and_protonation_composition(self, tmp_path):
        """Terminal prefix and protonation-state renaming compose correctly."""
        import MDAnalysis as mda

        pdb_path = self._make_terminal_protonation_pdb(tmp_path)
        u = mda.Universe(str(pdb_path))
        env_atoms = u.select_atoms("all")

        ens_parser = EnsembleParser()
        term_map = ens_parser._terminal_resname_map(u)
        prot_map = ens_parser._protonation_resname_map(env_atoms)

        # Chain A: N-terminal HIS (resindex 0) + HD1 → term_map should give NHIS
        # and prot_map should give HSD
        assert term_map[0] == "NHIS"
        assert prot_map[0] == "HSD"
        composed = ens_parser._prefixed_resname(prot_map[0], str(term_map[0])[0])
        assert composed == "NHSD"

        # Chain B: C-terminal HIS (should be last residue) + HE2
        # resindex 2 = ALA, resindex 3 = HIS (C-term via fallback OT1)
        assert term_map[3] == "CHIS"
        assert prot_map[3] == "HSE"
        composed = ens_parser._prefixed_resname(prot_map[3], str(term_map[3])[0])
        assert composed == "CHSE"

        # Chain C: N-terminal HIS (resindex 4) + both HD1/HE2 → NHSP
        assert term_map[4] == "NHIS"
        assert prot_map[4] == "HSP"
        composed = ens_parser._prefixed_resname(prot_map[4], str(term_map[4])[0])
        assert composed == "NHSP"

        # Chain D: N-terminal HIS (resindex 6) without HD1/HE2
        assert term_map[6] == "NHIS"
        assert 6 not in prot_map

        # Chain E: N-terminal GLU (resindex 8) + HE1/HE2 → NGLH
        assert term_map[8] == "NGLU"
        assert prot_map[8] == "GLH"
        composed = ens_parser._prefixed_resname(prot_map[8], str(term_map[8])[0])
        assert composed == "NGLH"

        # Chain F: C-terminal ASP (resindex 11) + HD2 → CASH
        assert term_map[11] == "CASP"
        assert prot_map[11] == "ASH"
        composed = ens_parser._prefixed_resname(prot_map[11], str(term_map[11])[0])
        assert composed == "CASH"

        # Chain G: N-terminal CYS (resindex 12) without HG → NCYX
        assert term_map[12] == "NCYS"
        assert prot_map[12] == "CYX"
        composed = ens_parser._prefixed_resname(prot_map[12], str(term_map[12])[0])
        assert composed == "NCYX"

        # The TIP3 water residue is present in the PDB but is neither terminal
        # protein nor a supported protonation-state residue.
        assert 14 not in term_map
        assert 14 not in prot_map

    def test_structures_applies_terminal_protonation_to_pe_resnames(self,
                                                                    tmp_path):
        """structures() emits composed terminal/protonation residue names."""

        pdb_path = self._make_terminal_protonation_pdb(tmp_path)
        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        snapshots = ens_parser.structures(
            pdb_file=pdb_path,
            qm_region="resname TIP3",
            pe_cutoff=150.0,
        )

        assert len(snapshots) == 1
        assert snapshots[0]["qm_atom_names"].tolist() == ["OH2", "H1", "H2"]
        pe_resnames = np.asarray(snapshots[0]["pe_resnames"], dtype=object)
        pe_resindices = np.asarray(snapshots[0]["pe_resindices"], dtype=int)
        resnames_by_resindex = {
            int(resindex): {
                str(resname)
                for resname in pe_resnames[pe_resindices == resindex]
            }
            for resindex in np.unique(pe_resindices)
        }

        assert resnames_by_resindex == {
            0: {"NHSD"},
            1: {"CALA"},
            2: {"NALA"},
            3: {"CHSE"},
            4: {"NHSP"},
            5: {"CALA"},
            6: {"NHIS"},
            7: {"CALA"},
            8: {"NGLH"},
            9: {"CALA"},
            10: {"NALA"},
            11: {"CASH"},
            12: {"NCYX"},
            13: {"CALA"},
        }
