"""
Pytest test suite for catcif_tools.

Biopython-dependent tests (TestGetSequence, TestChainTokenLengths,
TestConversion, TestCatcifsplitBio) are auto-skipped when Bio is not installed.
"""

import importlib.util
import os
import sys
import zlib

import pytest

_HAS_BIOPYTHON = importlib.util.find_spec('Bio') is not None


# ---------------------------------------------------------------------------
# Autouse fixture: clear the index cache between tests
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _clear_cache():
    from catcif_tools.cache import clear_cache
    clear_cache()
    yield
    clear_cache()


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _make_catcif(tmp_path, content, name='test.catcif'):
    p = tmp_path / name
    p.write_text(content)
    return p


def _data_lines(text):
    return [l.strip() for l in text.splitlines() if l.startswith('data_')]


# ---------------------------------------------------------------------------
# TestSettings
# ---------------------------------------------------------------------------

class TestSettings:

    def test_default_values(self):
        from catcif_tools.settings import catcif_settings
        assert catcif_settings.fast_cache_num == 3
        assert catcif_settings.fast_cache_time == 5000
        assert catcif_settings.slow_cache_num == 10
        assert catcif_settings.slow_cache_time == 100000
        assert catcif_settings.max_caches == 10

    def test_singleton(self):
        from catcif_tools.settings import catcif_settings
        import catcif_tools.settings as s2
        assert catcif_settings is s2.catcif_settings


# ---------------------------------------------------------------------------
# TestBuildIndex
# ---------------------------------------------------------------------------

class TestBuildIndex:

    def test_dedup_tag_first_occurrence(self):
        from catcif_tools.build_index import _dedup_tag
        seen = {}
        assert _dedup_tag('foo', seen) == 'foo'

    def test_dedup_tag_second_occurrence(self):
        from catcif_tools.build_index import _dedup_tag
        seen = {}
        _dedup_tag('foo', seen)
        assert _dedup_tag('foo', seen) == 'foo_1'

    def test_dedup_tag_third_occurrence(self):
        from catcif_tools.build_index import _dedup_tag
        seen = {}
        _dedup_tag('foo', seen)
        _dedup_tag('foo', seen)
        assert _dedup_tag('foo', seen) == 'foo_2'

    def test_add_runtime_fields(self):
        from catcif_tools.build_index import _add_runtime_fields
        index = {
            'version': 2,
            'orig_tags': ['a', 'b', 'c'],
            'index': {'a': {'o': 0}, 'b': {'o': 10}, 'c': {'o': 20}},
        }
        _add_runtime_fields(index)
        assert index['index']['a']['idx'] == 0
        assert index['index']['b']['idx'] == 1
        assert index['index']['c']['idx'] == 2

    def test_build_deduplicates_tags(self, tmp_path):
        from catcif_tools.build_index import build_catcif_index
        p = _make_catcif(tmp_path,
            'data_foo\n_cell.length_a 1.0\n'
            'data_bar\n_cell.length_a 2.0\n'
            'data_foo\n_cell.length_a 3.0\n')
        idx = build_catcif_index(str(p))
        assert list(idx['index'].keys()) == ['foo', 'bar', 'foo_1']
        assert idx['orig_tags'] == ['foo', 'bar', 'foo']

    def test_idx_absent_on_disk(self, tmp_path):
        import json
        from catcif_tools.build_index import build_catcif_index
        p = _make_catcif(tmp_path, 'data_foo\n_cell.length_a 1.0\n')
        build_catcif_index(str(p))
        raw = json.loads((tmp_path / 'test.catcif.idx').read_text())
        for entry in raw['index'].values():
            assert 'idx' not in entry

    def test_idx_present_in_memory(self, tmp_path):
        from catcif_tools.build_index import build_catcif_index
        p = _make_catcif(tmp_path,
            'data_a\n_cell.length_a 1.0\n'
            'data_b\n_cell.length_a 2.0\n')
        idx = build_catcif_index(str(p))
        for i, entry in enumerate(idx['index'].values()):
            assert entry['idx'] == i

    def test_no_tags_key_in_memory(self, tmp_path):
        from catcif_tools.build_index import build_catcif_index
        p = _make_catcif(tmp_path, 'data_foo\n_cell.length_a 1.0\n')
        assert 'tags' not in build_catcif_index(str(p))


# ---------------------------------------------------------------------------
# TestStructure
# ---------------------------------------------------------------------------

class TestStructure:

    def test_rename_header(self):
        from catcif_tools.structure import rename_structure
        result = rename_structure('data_old\n_entry.id old\n_cell.length_a 1.0\n', 'new')
        assert result.startswith('data_new\n')

    def test_rename_entry_id(self):
        from catcif_tools.structure import rename_structure
        result = rename_structure('data_old\n_entry.id old\n_cell.length_a 1.0\n', 'new')
        assert '_entry.id new\n' in result

    def test_rename_quoted_entry_id(self):
        from catcif_tools.structure import rename_structure
        result = rename_structure("data_foo\n_entry.id   'foo'\n_cell.length_a 1.0\n", 'bar')
        assert '_entry.id   bar' in result

    def test_rename_no_entry_id(self):
        from catcif_tools.structure import rename_structure
        result = rename_structure('data_foo\n_cell.length_a 1.0\n', 'bar')
        assert result.startswith('data_bar\n')

    def test_compress_gzip_magic(self):
        from catcif_tools.structure import compress_structure
        result = compress_structure('data_foo\n_entry.id foo\n')
        assert isinstance(result, bytes)
        assert result[:2] == b'\x1f\x8b'

    def test_compress_round_trip(self):
        from catcif_tools.structure import compress_structure
        s = 'data_foo\n_entry.id foo\n_cell.length_a 1.0\n'
        decompressed = zlib.decompress(compress_structure(s), wbits=47).decode('ascii')
        assert decompressed == s

    def test_get_all_structures_order(self, tmp_path):
        from catcif_tools.structure import get_all_structures
        p = _make_catcif(tmp_path,
            'data_a\n_entry.id a\n\ndata_b\n_entry.id b\n\ndata_c\n_entry.id c\n')
        tags = [tag for _, tag in get_all_structures(str(p))]
        assert tags == ['a', 'b', 'c']

    def test_get_all_structures_renames_duplicate(self, tmp_path):
        from catcif_tools.structure import get_all_structures
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_foo\n_entry.id foo\n')
        results = list(get_all_structures(str(p)))
        assert [tag for _, tag in results] == ['foo', 'foo_1']
        assert results[1][0].startswith('data_foo_1\n')

    def test_get_all_structures_preserve_tags(self, tmp_path):
        from catcif_tools.structure import get_all_structures
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_foo\n_entry.id foo\n')
        tags = [tag for _, tag in get_all_structures(str(p), preserve_tags=True)]
        assert tags == ['foo', 'foo']

    def test_get_all_structures_preserve_tags_text_unchanged(self, tmp_path):
        from catcif_tools.structure import get_all_structures
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_foo\n_entry.id foo\n')
        results = list(get_all_structures(str(p), preserve_tags=True))
        # Both structure texts should be unrenamed — data_foo, not data_foo_1
        assert results[0][0].startswith('data_foo\n')
        assert results[1][0].startswith('data_foo\n')

    def test_get_structure_by_tag(self, tmp_path):
        from catcif_tools.structure import get_structure
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_bar\n_entry.id bar\n')
        assert get_structure('bar', catcif_file=str(p)).startswith('data_bar\n')

    def test_get_structure_path_tag(self, tmp_path):
        from catcif_tools.structure import get_structure
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_bar\n_entry.id bar\n')
        assert get_structure(f'{p}:bar').startswith('data_bar\n')

    def test_get_structures_preserves_input_order(self, tmp_path):
        from catcif_tools.structure import get_structures
        p = _make_catcif(tmp_path,
            'data_a\n_entry.id a\n\ndata_b\n_entry.id b\n\ndata_c\n_entry.id c\n')
        # get_structures preserves the caller's requested order
        results = list(get_structures(['c', 'a'], catcif_file=str(p)))
        assert results[0].startswith('data_c\n')
        assert results[1].startswith('data_a\n')


# ---------------------------------------------------------------------------
# TestScores
# ---------------------------------------------------------------------------

class TestScores:

    BASE = 'data_foo\n_entry.id foo\n_cell.length_a 63.150\n'

    def test_write_scores_inserts_block(self):
        from catcif_tools.scores import write_scores
        result = write_scores(self.BASE, {'rmsd': 1.5, 'score': -10.0})
        assert '_catcif_scores.rmsd' in result
        assert '_catcif_scores.score' in result

    def test_get_scores_round_trip(self):
        from catcif_tools.scores import write_scores, get_scores
        s = write_scores(self.BASE, {'total_score': -123.45, 'rmsd': 1.234})
        assert get_scores(s) == {'total_score': '-123.45', 'rmsd': '1.234'}

    def test_write_scores_replaces_existing(self):
        from catcif_tools.scores import write_scores, get_scores
        s1 = write_scores(self.BASE, {'total_score': -1.0, 'rmsd': 1.0})
        s2 = write_scores(s1, {'total_score': 99.0, 'new_score': 'hello'})
        scores = get_scores(s2)
        assert scores == {'total_score': '99.0', 'new_score': 'hello'}
        assert 'rmsd' not in scores

    def test_write_scores_remove(self):
        from catcif_tools.scores import write_scores, get_scores
        s1 = write_scores(self.BASE, {'rmsd': 1.5})
        s2 = write_scores(s1, {})
        assert get_scores(s2) == {}

    def test_write_scores_noop_on_empty(self):
        from catcif_tools.scores import write_scores
        assert write_scores(self.BASE, {}) == self.BASE

    def test_parse_score_file(self, tmp_path):
        from catcif_tools.scores import write_scores, parse_score_file
        content = (
            write_scores('data_a\n_entry.id a\n', {'rmsd': 1.5, 'score': -10.0}) +
            write_scores('data_b\n_entry.id b\n', {'rmsd': 2.0, 'score': -5.0}) +
            write_scores('data_a\n_entry.id a\n', {'rmsd': 3.0, 'score': -1.0})
        )
        p = _make_catcif(tmp_path, content)
        results = list(parse_score_file(str(p)))
        assert results[0] == {'rmsd': '1.5', 'score': '-10.0', 'tag': 'a'}
        assert results[1] == {'rmsd': '2.0', 'score': '-5.0', 'tag': 'b'}
        assert results[2] == {'rmsd': '3.0', 'score': '-1.0', 'tag': 'a_1'}


# ---------------------------------------------------------------------------
# TestCatcifTools
# ---------------------------------------------------------------------------

class TestCatcifTools:

    def test_to_catcif_string_renames(self):
        from catcif_tools.catcif_tools import to_catcif_string
        s = to_catcif_string('data_old\n_entry.id old\n', 'new')
        assert s.startswith('data_new\n')
        assert '_entry.id new\n' in s

    def test_to_catcif_string_no_header_raises(self):
        from catcif_tools.catcif_tools import to_catcif_string
        with pytest.raises(ValueError, match='data_'):
            to_catcif_string('_cell.length_a 1.0\n', 'foo')

    def test_to_catcif_string_add_header(self):
        from catcif_tools.catcif_tools import to_catcif_string
        s = to_catcif_string('_cell.length_a 1.0\n', 'foo', add_header=True)
        assert s.startswith('data_foo\n')

    def test_to_catcif_string_with_scores(self):
        from catcif_tools.catcif_tools import to_catcif_string
        s = to_catcif_string('data_old\n_entry.id old\n', 'mytag',
                             scores={'rmsd': 1.5})
        assert '_catcif_scores.rmsd' in s

    def test_to_catcif_string_compress(self):
        from catcif_tools.catcif_tools import to_catcif_string
        result = to_catcif_string('data_old\n_entry.id old\n', 'mytag',
                                  compress=True)
        assert isinstance(result, bytes)
        assert result[:2] == b'\x1f\x8b'

    def test_append_to_catcif_file(self, tmp_path):
        from catcif_tools.catcif_tools import append_to_catcif_file
        from catcif_tools.structure import get_all_structures
        p = tmp_path / 'out.catcif'
        append_to_catcif_file(str(p), 'data_a\n_entry.id a\n', 'a')
        append_to_catcif_file(str(p), 'data_b\n_entry.id b\n', 'b')
        tags = [tag for _, tag in get_all_structures(str(p))]
        assert tags == ['a', 'b']

    def test_append_mixed_compressed_plain_round_trip(self, tmp_path):
        from catcif_tools.catcif_tools import append_to_catcif_file
        from catcif_tools.structure import get_all_structures
        p = tmp_path / 'mixed.catcif'
        append_to_catcif_file(str(p), 'data_a\n_entry.id a\n_cell.length_a 1.0\n', 'a')
        append_to_catcif_file(str(p), 'data_b\n_entry.id b\n_cell.length_a 2.0\n', 'b',
                              compress=True)
        append_to_catcif_file(str(p), 'data_c\n_entry.id c\n_cell.length_a 3.0\n', 'c')
        results = list(get_all_structures(str(p)))
        assert [tag for _, tag in results] == ['a', 'b', 'c']
        assert results[1][0].startswith('data_b\n')

    def test_get_tags(self, tmp_path):
        from catcif_tools.catcif_tools import get_tags, get_tags_from_index
        from catcif_tools.cache import get_catcif_index
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_bar\n_entry.id bar\n')
        index = get_catcif_index(str(p))
        assert get_tags_from_index(index) == ['foo', 'bar']
        assert get_tags(str(p)) == ['foo', 'bar']


# ---------------------------------------------------------------------------
# TestPath
# ---------------------------------------------------------------------------

class TestPath:

    def test_split_bare_tag(self):
        from catcif_tools.path import split_catcif_tag
        path, tag = split_catcif_tag('mytag')
        assert path is None
        assert tag == 'mytag'

    def test_split_catcif_path_tag(self):
        from catcif_tools.path import split_catcif_tag
        path, tag = split_catcif_tag('/some/file.catcif:mytag')
        assert path == '/some/file.catcif'
        assert tag == 'mytag'

    def test_split_cif_path_tag(self):
        from catcif_tools.path import split_catcif_tag
        path, tag = split_catcif_tag('/some/file.cif:mytag')
        assert path == '/some/file.cif'
        assert tag == 'mytag'

    def test_split_preserves_colons_in_tag(self):
        from catcif_tools.path import split_catcif_tag
        path, tag = split_catcif_tag('/some/file.catcif:a:b:c')
        assert path == '/some/file.catcif'
        assert tag == 'a:b:c'

    def test_is_path_tag_true(self):
        from catcif_tools.path import is_catcif_path_tag
        assert is_catcif_path_tag('/some/file.catcif:mytag')
        assert is_catcif_path_tag('/some/file.cif:mytag')

    def test_is_path_tag_false(self):
        from catcif_tools.path import is_catcif_path_tag
        assert not is_catcif_path_tag('mytag')
        assert not is_catcif_path_tag('not_a_path')


# ---------------------------------------------------------------------------
# TestCache
# ---------------------------------------------------------------------------

class TestCache:

    def test_singleton(self):
        from catcif_tools.cache import catcif_cache
        import catcif_tools.cache as c2
        assert catcif_cache is c2.catcif_cache

    def test_get_catcif_index_returns_index(self, tmp_path):
        from catcif_tools.cache import get_catcif_index
        p = _make_catcif(tmp_path, 'data_foo\n_entry.id foo\n')
        index = get_catcif_index(str(p))
        assert 'foo' in index['index']
        assert index['version'] == 2

    def test_clear_cache(self, tmp_path):
        from catcif_tools.cache import get_catcif_index, clear_cache, catcif_cache
        p = _make_catcif(tmp_path, 'data_foo\n_entry.id foo\n')
        get_catcif_index(str(p))
        assert len(catcif_cache) > 0
        clear_cache()
        assert len(catcif_cache) == 0


# ---------------------------------------------------------------------------
# TestCatcifls
# ---------------------------------------------------------------------------

class TestCatcifls:

    def test_basic_list(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifls import main
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_bar\n_entry.id bar\n'
            '\ndata_foo\n_entry.id foo\n')
        monkeypatch.setattr(sys, 'argv', ['catcifls', str(p)])
        main()
        assert capsys.readouterr().out == 'foo\nbar\nfoo_1\n'

    def test_g_flag_original_names(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifls import main
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_bar\n_entry.id bar\n'
            '\ndata_foo\n_entry.id foo\n')
        monkeypatch.setattr(sys, 'argv', ['catcifls', '-g', str(p)])
        main()
        assert capsys.readouterr().out == 'foo\nbar\nfoo\n'

    def test_l_flag_prefixes_path(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifls import main
        p = _make_catcif(tmp_path, 'data_foo\n_entry.id foo\n')
        monkeypatch.setattr(sys, 'argv', ['catcifls', '-l', str(p)])
        main()
        out = capsys.readouterr().out
        assert out == f'{p.resolve()}:foo\n'

    def test_l_and_g_flags_combined(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifls import main
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n\ndata_bar\n_entry.id bar\n'
            '\ndata_foo\n_entry.id foo\n')
        monkeypatch.setattr(sys, 'argv', ['catcifls', '-l', '-g', str(p)])
        main()
        rp = str(p.resolve())
        assert capsys.readouterr().out.splitlines() == [
            f'{rp}:foo', f'{rp}:bar', f'{rp}:foo']


# ---------------------------------------------------------------------------
# TestCatcifslice
# ---------------------------------------------------------------------------

class TestCatcifslice:

    def test_slice_single_tag(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifslice import main
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n_cell.length_a 1.0\n'
            'data_bar\n_entry.id bar\n_cell.length_a 2.0\n')
        monkeypatch.setattr(sys, 'argv', ['catcifslice', str(p), 'foo'])
        main()
        out = capsys.readouterr().out
        assert 'data_foo' in out
        assert 'data_bar' not in out

    def test_missing_tag_warns_on_stderr(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifslice import main
        p = _make_catcif(tmp_path, 'data_foo\n_entry.id foo\n')
        monkeypatch.setattr(sys, 'argv', ['catcifslice', str(p), 'missing'])
        with pytest.raises(SystemExit):
            main()
        assert 'Unable to find' in capsys.readouterr().err

    def test_no_tags_exits(self, tmp_path, monkeypatch):
        from catcif_tools.catcifslice import main
        p = _make_catcif(tmp_path, 'data_foo\n_entry.id foo\n')
        monkeypatch.setattr(sys, 'argv', ['catcifslice', str(p)])
        with pytest.raises(SystemExit):
            main()


# ---------------------------------------------------------------------------
# TestCatcifextract
# ---------------------------------------------------------------------------

class TestCatcifextract:

    def test_extract_writes_cif_files(self, tmp_path, monkeypatch):
        from catcif_tools.catcifextract import main
        p = _make_catcif(tmp_path,
            'data_foo\n_entry.id foo\n_cell.length_a 1.0\n'
            'data_bar\n_entry.id bar\n_cell.length_a 2.0\n')
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv',
                            ['catcifextract', '-o', str(out), str(p), 'foo', 'bar'])
        main()
        assert (out / 'foo.cif').read_text().startswith('data_foo')
        assert (out / 'bar.cif').read_text().startswith('data_bar')

    def test_extract_z_writes_gz_files(self, tmp_path, monkeypatch):
        from catcif_tools.catcifextract import main
        p = _make_catcif(tmp_path, 'data_foo\n_entry.id foo\n_cell.length_a 1.0\n')
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv',
                            ['catcifextract', '-z', '-o', str(out), str(p), 'foo'])
        main()
        data = (out / 'foo.cif.gz').read_bytes()
        assert data[:2] == b'\x1f\x8b'


# ---------------------------------------------------------------------------
# TestCatciffromfiles
# ---------------------------------------------------------------------------

class TestCatciffromfiles:

    def test_tag_from_path(self):
        from catcif_tools.catciffromfiles import _tag_from_path
        assert _tag_from_path('foo/bar/myprotein.cif') == 'myprotein'
        assert _tag_from_path('foo/bar/myprotein.cif.gz') == 'myprotein'

    def test_from_cif_file(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catciffromfiles import main
        cif = tmp_path / 'myprotein.cif'
        cif.write_text('data_myprotein\n_entry.id myprotein\n_cell.length_a 1.0\n')
        monkeypatch.setattr(sys, 'argv', ['catciffromfiles', str(cif)])
        main()
        assert 'data_myprotein' in capsys.readouterr().out

    def test_unsupported_extension_exits(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catciffromfiles import main
        bad = tmp_path / 'file.txt'
        bad.write_text('not a cif')
        monkeypatch.setattr(sys, 'argv', ['catciffromfiles', str(bad)])
        with pytest.raises(SystemExit):
            main()


# ---------------------------------------------------------------------------
# TestCatcifscorefile
# ---------------------------------------------------------------------------

class TestCatcifscorefile:

    def test_writes_sc_file(self, tmp_path, monkeypatch):
        from catcif_tools.catcifscorefile import main
        from catcif_tools.scores import write_scores
        content = (
            write_scores('data_foo\n_entry.id foo\n', {'rmsd': 1.5, 'score': -10.0}) +
            write_scores('data_bar\n_entry.id bar\n', {'rmsd': 2.0, 'score': -5.0})
        )
        p = _make_catcif(tmp_path, content)
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv', ['catcifscorefile', str(p), '-o', str(out)])
        main()
        lines = (out / 'test.sc').read_text().splitlines()
        assert lines[0].split() == ['rmsd', 'score', 'tag']
        assert lines[1].split() == ['1.5', '-10.0', 'foo']
        assert lines[2].split() == ['2.0', '-5.0', 'bar']

    def test_scorefile_deduplicates_tags(self, tmp_path, monkeypatch):
        from catcif_tools.catcifscorefile import main
        from catcif_tools.scores import write_scores
        content = (
            write_scores('data_foo\n_entry.id foo\n', {'rmsd': 1.5, 'score': -10.0}) +
            write_scores('data_bar\n_entry.id bar\n', {'rmsd': 2.0, 'score': -5.0}) +
            write_scores('data_foo\n_entry.id foo\n', {'rmsd': 3.0, 'score': -1.0})
        )
        p = _make_catcif(tmp_path, content)
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv', ['catcifscorefile', str(p), '-o', str(out)])
        main()
        lines = (out / 'test.sc').read_text().splitlines()
        assert lines[0].split() == ['rmsd', 'score', 'tag']
        assert lines[1].split() == ['1.5', '-10.0', 'foo']
        assert lines[2].split() == ['2.0', '-5.0', 'bar']
        assert lines[3].split() == ['3.0', '-1.0', 'foo_1']

    def test_writes_csv_with_c_flag(self, tmp_path, monkeypatch):
        from catcif_tools.catcifscorefile import main
        from catcif_tools.scores import write_scores
        content = write_scores('data_foo\n_entry.id foo\n', {'rmsd': 1.5})
        p = _make_catcif(tmp_path, content)
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv',
                            ['catcifscorefile', str(p), '-o', str(out), '-c'])
        main()
        lines = (out / 'test.csv').read_text().splitlines()
        assert lines[0] == 'rmsd,tag'


# ---------------------------------------------------------------------------
# TestCatcifrename
# ---------------------------------------------------------------------------

class TestCatcifrename:

    def test_rename_basic(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifrename import main
        p = _make_catcif(tmp_path,
            'data_a\n_entry.id a\n_cell.length_a 1.0\n'
            'data_b\n_entry.id b\n_cell.length_a 2.0\n')
        monkeypatch.setattr(sys, 'argv', ['catcifrename', str(p), 'alpha', 'beta'])
        main()
        out = capsys.readouterr().out
        assert 'data_alpha' in out
        assert 'data_beta' in out
        assert 'data_a\n' not in out

    def test_wrong_count_exits(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifrename import main
        p = _make_catcif(tmp_path,
            'data_a\n_entry.id a\n\ndata_b\n_entry.id b\n')
        monkeypatch.setattr(sys, 'argv', ['catcifrename', str(p), 'only_one'])
        with pytest.raises(SystemExit):
            main()
        assert 'does not match' in capsys.readouterr().err

    def test_broadcast_mode(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifrename import main
        p = _make_catcif(tmp_path, 'data_tmpl\n_entry.id tmpl\n_cell.length_a 1.0\n')
        monkeypatch.setattr(sys, 'argv',
                            ['catcifrename', '-b', str(p), 'x1', 'x2', 'x3'])
        main()
        out = capsys.readouterr().out
        assert out.count('data_x') == 3


# ---------------------------------------------------------------------------
# TestCatcifsequence
# ---------------------------------------------------------------------------

class TestCatcifsequence:

    _CIF = (
        'data_foo\n#\nloop_\n'
        '_atom_site.label_entity_id\n'
        '_atom_site.label_asym_id\n'
        '_atom_site.label_comp_id\n'
        '_atom_site.label_seq_id\n'
        '1 A ALA 1\n1 A GLY 2\n1 A TRP 3\n'
        '2 B PHE 1\n2 B VAL 2\n'
    )

    def test_sequence_output(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifsequence import main
        p = _make_catcif(tmp_path, self._CIF)
        monkeypatch.setattr(sys, 'argv', ['catcifsequence', str(p)])
        main()
        assert capsys.readouterr().out.strip() == 'AGW FV foo'

    def test_c_flag_prefixes_chain_id(self, tmp_path, monkeypatch, capsys):
        from catcif_tools.catcifsequence import main
        p = _make_catcif(tmp_path, self._CIF)
        monkeypatch.setattr(sys, 'argv', ['catcifsequence', '-c', str(p)])
        main()
        assert capsys.readouterr().out.strip() == 'A:AGW B:FV foo'


# ---------------------------------------------------------------------------
# TestCatcifsplit
# ---------------------------------------------------------------------------

class TestCatcifsplit:

    def test_suffix_generator_start(self):
        from catcif_tools.catcifsplit import _suffix_gen
        gen = _suffix_gen()
        assert [next(gen) for _ in range(3)] == ['xaa', 'xab', 'xac']

    def test_suffix_generator_transition(self):
        from catcif_tools.catcifsplit import _suffix_gen
        gen = _suffix_gen()
        for _ in range(675):
            next(gen)
        assert next(gen) == 'xzz'
        assert next(gen) == 'xaaa'

    def test_plain_split(self, tmp_path, monkeypatch):
        from catcif_tools.catcifsplit import main
        content = ''.join(
            f'data_s{i}\n_entry.id s{i}\n' for i in range(5))
        p = _make_catcif(tmp_path, content)
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv', ['catcifsplit', str(p), '2', '-o', str(out)])
        main()
        files = sorted(os.listdir(str(out)))
        assert files == ['xaa.catcif', 'xab.catcif', 'xac.catcif']
        assert _data_lines((out / 'xaa.catcif').read_text()) == ['data_s0', 'data_s1']
        assert _data_lines((out / 'xac.catcif').read_text()) == ['data_s4']

    def test_shuffle_split(self, tmp_path, monkeypatch):
        from catcif_tools.catcifsplit import main
        content = ''.join(
            f'data_s{i}\n_entry.id s{i}\n' for i in range(6))
        p = _make_catcif(tmp_path, content)
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv',
                            ['catcifsplit', '-s', str(p), '2', '-o', str(out)])
        main()
        assert sorted(os.listdir(str(out))) == ['xaa.catcif', 'xab.catcif', 'xac.catcif']
        assert _data_lines((out / 'xaa.catcif').read_text()) == ['data_s0', 'data_s3']
        assert _data_lines((out / 'xab.catcif').read_text()) == ['data_s1', 'data_s4']
        assert _data_lines((out / 'xac.catcif').read_text()) == ['data_s2', 'data_s5']


# ---------------------------------------------------------------------------
# Biopython-dependent tests
# ---------------------------------------------------------------------------

_BIO_CIF = (
    'data_test\n#\nloop_\n'
    '_atom_site.label_entity_id\n'
    '_atom_site.label_asym_id\n'
    '_atom_site.label_comp_id\n'
    '_atom_site.label_seq_id\n'
    '1 A ALA 1\n'
    '1 A ALA 1\n'   # duplicate atom — same residue
    '1 A GLY 2\n'
    '1 A TRP 5\n'   # gap in seq_id — NOT a new chain
    '1 A MSE 6\n'   # non-standard
    '2 B PHE 1\n'
    '2 B VAL 2\n'
    '3 A HIS 1\n'   # same chain_id A but new entity_id — IS a new chain
)

_TOKEN_CIF = (
    'data_test\n#\nloop_\n'
    '_atom_site.label_entity_id\n'
    '_atom_site.label_asym_id\n'
    '_atom_site.label_comp_id\n'
    '_atom_site.label_seq_id\n'
    '_atom_site.auth_seq_id\n'
    '_atom_site.type_symbol\n'
    '1 A ALA 1 1 N\n'
    '1 A ALA 1 1 C\n'
    '1 A ALA 1 1 H\n'   # hydrogen — excluded from non-canonical count
    '1 A GLY 2 2 N\n'
    '1 A MSE 3 3 C\n'
    '1 A MSE 3 3 SE\n'
    '1 A MSE 3 3 H\n'   # hydrogen
    '2 B LIG . 101 C\n'
    '2 B LIG . 101 N\n'
    '2 B LIG . 101 H\n'
    '2 B LIG . 102 C\n'
)


class TestGetSequence:

    @pytest.fixture(autouse=True)
    def _require_bio(self):
        pytest.importorskip('Bio')

    def test_three_chains_returned(self):
        from catcif_tools.biology import get_sequence
        assert len(get_sequence(_BIO_CIF)) == 3

    def test_canonical_sequence(self):
        from catcif_tools.biology import get_sequence
        result = get_sequence(_BIO_CIF)
        assert result[0]['sequence'] == 'AGWX'

    def test_duplicate_atom_deduped(self):
        from catcif_tools.biology import get_sequence
        # ALA appears twice at seq_id 1; should count once
        result = get_sequence(_BIO_CIF)
        assert len(result[0]['sequence']) == 4  # not 5

    def test_gap_does_not_split_chain(self):
        from catcif_tools.biology import get_sequence
        # seq_id jumps from 2 to 5 — still one chain entry
        result = get_sequence(_BIO_CIF)
        assert result[0]['chain_id'] == 'A'
        assert result[0]['entity_id'] == '1'

    def test_nonstandard_residue_is_x(self):
        from catcif_tools.biology import get_sequence
        result = get_sequence(_BIO_CIF)
        assert result[0]['sequence'][3] == 'X'

    def test_entity_change_starts_new_chain(self):
        from catcif_tools.biology import get_sequence
        result = get_sequence(_BIO_CIF)
        assert result[2] == {'chain_id': 'A', 'entity_id': '3', 'sequence': 'H'}

    def test_second_chain(self):
        from catcif_tools.biology import get_sequence
        result = get_sequence(_BIO_CIF)
        assert result[1] == {'chain_id': 'B', 'entity_id': '2', 'sequence': 'FV'}


class TestChainTokenLengths:

    @pytest.fixture(autouse=True)
    def _require_bio(self):
        pytest.importorskip('Bio')

    def test_two_chains(self):
        from catcif_tools.biology import chain_token_lengths
        result = chain_token_lengths(_TOKEN_CIF)
        assert len(result) == 2

    def test_canonical_aa_is_one_token(self):
        from catcif_tools.biology import chain_token_lengths
        result = chain_token_lengths(_TOKEN_CIF)
        # ALA=1, GLY=1, MSE=2 non-H (C + SE) → 4
        assert result[0] == {'chain_id': 'A', 'entity_id': '1', 'tokens': 4}

    def test_hydrogen_excluded(self):
        from catcif_tools.biology import chain_token_lengths
        result = chain_token_lengths(_TOKEN_CIF)
        # MSE has 3 atoms but 1 H → 2 non-H tokens, not 3
        assert result[0]['tokens'] == 4  # would be 5 if H were counted

    def test_ligand_groups_by_auth_seq_id(self):
        from catcif_tools.biology import chain_token_lengths
        result = chain_token_lengths(_TOKEN_CIF)
        # LIG@101: C+N (H excluded) = 2; LIG@102: C = 1 → 3
        assert result[1] == {'chain_id': 'B', 'entity_id': '2', 'tokens': 3}


class TestConversion:

    _PDB = (
        'ATOM      1  CA  ALA A   1       1.000   2.000   3.000'
        '  1.00  0.00           C\nEND\n'
    )

    @pytest.fixture(autouse=True)
    def _require_bio(self):
        pytest.importorskip('Bio')

    def test_returns_string(self):
        from catcif_tools.conversion import pdb_to_cif
        assert isinstance(pdb_to_cif(self._PDB, 'myprotein'), str)

    def test_has_data_header(self):
        from catcif_tools.conversion import pdb_to_cif
        result = pdb_to_cif(self._PDB, 'myprotein')
        assert result.splitlines()[0] == 'data_myprotein'

    def test_process_pdb_file_round_trip(self, tmp_path):
        from catcif_tools.catciffromfiles import _process_file
        from catcif_tools.structure import get_all_structures
        pdb_file = tmp_path / 'myprotein.pdb'
        pdb_file.write_text(
            'ATOM      1  CA  ALA A   1       1.000   2.000   3.000'
            '  1.00  0.00           C\nEND\n'
        )
        result = _process_file(str(pdb_file))
        assert result.startswith('data_myprotein')
        catcif = tmp_path / 'out.catcif'
        catcif.write_text(result)
        tags = [tag for _, tag in get_all_structures(str(catcif))]
        assert tags == ['myprotein']

    def test_import_error_without_biopython(self, monkeypatch):
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name.startswith('Bio'):
                raise ImportError('mocked missing biopython')
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, '__import__', mock_import)
        from catcif_tools.conversion import pdb_to_cif_bio
        with pytest.raises(ImportError, match='biopython'):
            pdb_to_cif_bio(self._PDB, 'test')


class TestCatcifsplitBio:
    """Tests for -t and -b modes which require biopython."""

    @pytest.fixture(autouse=True)
    def _require_bio(self):
        pytest.importorskip('Bio')

    _CIF_A = (
        'data_TAG\n#\nloop_\n'
        '_atom_site.label_entity_id\n_atom_site.label_asym_id\n'
        '_atom_site.label_comp_id\n_atom_site.label_seq_id\n'
        '_atom_site.auth_seq_id\n_atom_site.type_symbol\n'
        '1 A ALA 1 1 N\n1 A ALA 2 2 N\n1 A ALA 3 3 N\n'
    )
    _CIF_AB = (
        'data_TAG\n#\nloop_\n'
        '_atom_site.label_entity_id\n_atom_site.label_asym_id\n'
        '_atom_site.label_comp_id\n_atom_site.label_seq_id\n'
        '_atom_site.auth_seq_id\n_atom_site.type_symbol\n'
        '1 A ALA 1 1 N\n1 A ALA 2 2 N\n'
        '2 B GLY 1 1 N\n2 B GLY 2 2 N\n'
    )

    def _make_multi_catcif(self, tmp_path):
        structs = [
            self._CIF_A.replace('data_TAG', f'data_s{i}')
            for i in range(3)
        ] + [
            self._CIF_AB.replace('data_TAG', f'data_s{i}')
            for i in range(3, 7)
        ]
        return _make_catcif(tmp_path, ''.join(structs))

    def test_t_mode_groups_by_chain_tokens(self, tmp_path, monkeypatch):
        from catcif_tools.catcifsplit import main
        p = self._make_multi_catcif(tmp_path)
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv',
                            ['catcifsplit', '-t', str(p), '2', '-o', str(out)])
        main()
        files = sorted(os.listdir(str(out)))
        assert any(f.startswith('len3') for f in files)
        assert any(f.startswith('len2_2') for f in files)

    def test_b_mode_groups_by_bucketed_total(self, tmp_path, monkeypatch):
        from catcif_tools.catcifsplit import main
        p = self._make_multi_catcif(tmp_path)
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv',
                            ['catcifsplit', '-b', '--bucket-step', '2',
                             str(p), '2', '-o', str(out)])
        main()
        files = sorted(os.listdir(str(out)))
        # 3 tokens → bucketed to 2; 4 tokens stays 4
        assert any(f.startswith('len2') for f in files)
        assert any(f.startswith('len4') for f in files)
