// SPDX-License-Identifier: MIT — option form for the CRISPRuno settings (help texts from CRISPRunoCore.py argparse, CRISPRuno 0.1.22)
'use strict';
/* The settings textarea stays the source of truth: the form reads its values from CuSettings.parse() and every
 * change in the form rewrites the matching "key<TAB>value" line (a line is removed when the value equals the default). */
(function (root) {
const GROUPS = [
	['Input files', 'File names only; the files themselves are chosen under "2. Input files".', [
		['fastq_r1', 'file', 'Input FASTQ R1. Reads in this file are primed from the primer sequence.'],
		['fastq_r2', 'file', 'Input FASTQ R2 (optional).'],
		['fastq_umi', 'file', 'Input FASTQ with the UMI sequences (optional).'],
		['cut_region_annotation_file', 'file', 'BED file (chr, start, end, region_name) used to annotate cut sites that are not otherwise annotated (optional).'],
	]],
	['Primer and trimming', '', [
		['primer_seq', 'seq', 'Sequence of the primer (required).'],
		['guide_sequences', 'list', 'Spacer sequences of the guides without the PAM, oriented so the PAM would follow the spacer. Separate several guides with spaces.'],
		['primer_in_r2', 'bool', 'The primer is in R2. By default the primer is required to be in R1.'],
		['min_primer_aln_score', 'int', 'Minimum primer/origin alignment score for trimming.'],
		['min_primer_length', 'int', 'Minimum length of sequence required to match between the primer/origin and the read.'],
		['allow_indels_in_origin_aln', 'bool', 'Allow indels in the primer/origin alignment. By default a read with an indel in the primer/origin is discarded.'],
		['min_read_length', 'int', 'Minimum length of a read after all filtering.'],
		['transposase_adapter_seq', 'seq', 'Transposase adapter sequence trimmed from the reads (cutadapt -a).'],
	]],
	['Alignment filtering', 'minibwa MAPQ goes up to 60 (bowtie2: 42), so MAPQ thresholds are not on exactly the same scale as in the original.', [
		['min_alignment_quality_score', 'int', 'Minimum alignment quality (MAPQ) to accept an alignment.'],
		['arm_min_matched_start_bases', 'int', 'Number of bases at each end of the read alignment that must match exactly (no indels or mismatches). Reads failing this are "poor alignments".'],
		['arm_max_clipped_bases', 'int', 'Maximum number of clipped bases at the start of the alignment. Clipped bases are expected on one side only (e.g. from insertions); reads clipped more than this on both sides are poor alignments.'],
		['ignore_n', 'bool', 'Ignore N bases. By default N counts as a mismatch for arm_min_matched_start_bases.'],
		['suppress_poor_alignment_filter', 'bool', 'Keep reads with poor alignments in the final analysis and counts (excluded by default).'],
		['analyze_multimap_assignments', 'bool', 'Analyze reads that map to several locations and assign them to locations near known cut sites. Otherwise only the primary alignment is used.'],
		['max_multimap_alignments', 'int', 'Maximum number of multimapping alignments reported per read (with analyze_multimap_assignments).'],
	]],
	['R1/R2 support', '', [
		['r1_r2_support_max_distance', 'int', 'Maximum distance between R1 and R2 for the pair to count as supported by R2.'],
		['min_r2_alignment_quality_score', 'int', 'Minimum MAPQ of R2 for the pair to count as supported by R2.'],
		['suppress_r2_support_filter', 'bool', 'Keep reads without R2 support in the final analysis and counts (excluded by default).'],
	]],
	['UMI and deduplication', '', [
		['umi_regex', 'str', 'Pattern the UMI must match (IUPAC, e.g. NNWNNWNNN). Empty = no check.'],
		['dedup_input_on_UMI', 'bool', 'Deduplicate input reads on UMI before alignment.'],
		['suppress_dedup_on_aln_pos_and_UMI_filter', 'bool', 'Keep reads called duplicates by alignment position and UMI in the final analysis and counts.'],
		['dedup_by_final_cut_assignment_and_UMI', 'bool', 'Deduplicate on final cut assignment and UMI, so reads with the same UMI but different alignment start/stop are merged when assigned to the same cut.'],
		['min_umi_seen_to_keep_read', 'int', 'Minimum times a UMI/read combination must be seen to keep it.'],
		['write_UMI_counts', 'bool', 'Write a file with each UMI and its number of reads.'],
	]],
	['Cut sites', '', [
		['cleavage_offset', 'int', 'Cleavage position relative to the end of the spacer (Cas9: -3).'],
		['novel_cut_merge_distance', 'int', 'Novel cut sites within this distance (bp) of each other are merged into one site.'],
		['known_cut_merge_distance', 'int', 'Novel cut sites within this distance (bp) of a known/given/homologous site (not the origin) are merged into that site.'],
		['origin_cut_merge_distance', 'int', 'Reads aligned within this distance (bp) of the origin are merged into the origin.'],
		['suppress_homology_detection', 'bool', 'Skip the search for guide homology at novel cut sites.'],
		['guide_homology_max_gaps', 'int', 'Maximum gaps in the guide alignment for a site to count as homologous.'],
		['guide_homology_max_mismatches', 'int', 'Maximum mismatches in the guide alignment for a site to count as homologous.'],
		['cut_classification_annotations', 'annos', 'Custom labels for cut products, e.g. chr1:234:left:Custom_label. Separate several with spaces.'],
	]],
	['CRISPResso', '', [
		['run_crispresso_on_novel_sites', 'bool', 'Also run CRISPResso on novel cut sites (otherwise only on the origin and given on/off-targets).'],
		['crispresso_min_count', 'int', 'Minimum number of reads at a site for it to be analyzed by CRISPResso.'],
		['crispresso_max_indel_size', 'int', 'Maximum indel length (from the genome alignment) for a read to be analyzed by CRISPResso.'],
		['crispresso_min_aln_score', 'int', 'Minimum alignment score to the reference for CRISPResso quantification.'],
		['crispresso_quant_window_size', 'int', 'Number of bp on each side of the cut considered for edits.'],
	]],
	['Output', '', [
		['keep_intermediate', 'bool', 'Keep intermediate files.'],
		['write_discarded_read_info', 'bool', 'Write a file with information on discarded reads.'],
		['suppress_plots', 'bool', 'Do not make plots.'],
		['debug', 'bool', 'Print debug output.'],
	]],
];

const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const lineKey = (line) => (line.startsWith('#') ? null : line.split('#')[0].split('\t')[0].trim());

/* value of a setting as written in the settings file ('' = default / line removed) */
function toText(type, v, def) {
	if (type === 'bool') return v === def ? '' : (v ? 'True' : 'False');
	if (type === 'list' || type === 'annos') { const a = (v || []).map((x) => (type === 'annos' ? String(x).replace(/ /g, '__') : x)); return a.length ? a.join(' ') : ''; }
	if (v === null || v === undefined) return '';
	return String(v) === String(def) ? '' : String(v);
}
/* set / replace / remove the line of `key` in the settings text */
function setLine(text, key, value) {
	const lines = text.split('\n');
	const idx = lines.map((l, i) => (lineKey(l.replace(/\r$/, '')) === key ? i : -1)).filter((i) => i >= 0);
	if (value === '') return lines.filter((_, i) => !idx.includes(i)).join('\n');
	const nl = `${key}\t${value}`;
	if (idx.length) { lines[idx[0]] = nl; return lines.filter((_, i) => !idx.slice(1).includes(i)).join('\n'); }
	while (lines.length && lines[lines.length - 1].trim() === '') lines.pop();
	lines.push(nl);
	return lines.join('\n') + '\n';
}

/* build(container, { getText, setText, defaults }) -> { sync(settings) } */
function build(el, { getText, setText, defaults }) {
	const DEF = Object.fromEntries(defaults.map(([k, d, t]) => [k, { d, t }]));
	const fields = [];
	let html = '';
	for (const [g, note, opts] of GROUPS) {
		html += `<details class="optgrp"${g === 'Primer and trimming' || g === 'Alignment filtering' ? ' open' : ''}><summary>${esc(g)} <span class="badge" data-grp="${esc(g)}"></span></summary>` +
			(note ? `<p class="note">${esc(note)}</p>` : '') + '<div class="optgrid">';
		for (const [k, type, help] of opts) {
			const d = DEF[k] ? DEF[k].d : null, id = 'opt-' + k;
			fields.push({ k, type, d, g });
			const dtxt = type === 'bool' ? (d ? 'True' : 'False') : Array.isArray(d) ? (d.length ? d.join(' ') : 'none') : d === null ? 'none' : String(d);
			const input = type === 'bool' ? `<input type="checkbox" id="${id}">`
				: type === 'int' ? `<input type="number" id="${id}" step="1" placeholder="${esc(dtxt)}">`
				: `<input type="text" id="${id}" class="mono" spellcheck="false" placeholder="${esc(dtxt)}"${type === 'file' ? ' list="opt-files"' : ''}>`;
			html += `<label for="${id}" class="mono" title="${esc(k)}">${esc(k)}</label><div class="optin">${input} <button type="button" class="reset hidden" data-k="${esc(k)}" title="Reset to the default">↺</button>` +
				`<div class="note">${esc(help)} <span class="def">Default: <span class="mono">${esc(dtxt)}</span></span></div></div>`;
		}
		html += '</div></details>';
	}
	el.innerHTML = html + '<datalist id="opt-files"></datalist><div class="row"><button type="button" id="opt-reset-all">Reset all options to defaults</button><span class="note">Input files and primer/guide sequences are kept.</span></div>';

	const $f = (k) => el.querySelector('#opt-' + k);
	const write = (f) => {
		const inp = $f(f.k);
		let v;
		if (f.type === 'bool') v = inp.checked === !!f.d ? '' : (inp.checked ? 'True' : 'False');
		else { v = inp.value.trim(); if (f.type === 'int' && v !== '' && String(parseInt(v, 10)) === String(f.d)) v = ''; }
		if (f.type === 'annos') v = v.split(/\s+/).filter(Boolean).join(' ');
		setText(setLine(getText(), f.k, v));
	};
	for (const f of fields) {
		const inp = $f(f.k);
		inp.addEventListener(f.type === 'bool' ? 'change' : 'input', () => write(f));
	}
	el.addEventListener('click', (e) => {
		const b = e.target.closest('button.reset');
		if (b) { setText(setLine(getText(), b.dataset.k, '')); return; }
		if (e.target.id === 'opt-reset-all') {
			let t = getText();
			for (const f of fields) if (!['file', 'seq', 'list'].includes(f.type) || f.k === 'transposase_adapter_seq') t = setLine(t, f.k, '');
			setText(t);
		}
	});

	/* show the parsed settings in the form (skips the field being edited) */
	function sync(S, fileNames = []) {
		el.querySelector('#opt-files').innerHTML = fileNames.map((n) => `<option value="${esc(n)}">`).join('');
		const changed = {};
		for (const f of fields) {
			const inp = $f(f.k), v = S ? S[f.k] : f.d;
			const isDef = JSON.stringify(v) === JSON.stringify(f.d) || (v === null && f.d === null);
			if (document.activeElement !== inp) {
				if (f.type === 'bool') inp.checked = !!v;
				else if (f.type === 'list' || f.type === 'annos') inp.value = (v || []).map((x) => (f.type === 'annos' ? String(x).replace(/ /g, '__') : x)).join(' ');
				else inp.value = v === null || v === undefined ? '' : (f.type === 'int' && isDef ? '' : String(v));
			}
			const show = !isDef && !['file', 'seq', 'list'].includes(f.type) || (f.k === 'transposase_adapter_seq' && !isDef);
			el.querySelector(`button.reset[data-k="${f.k}"]`).classList.toggle('hidden', !show);
			inp.closest('.optin').previousElementSibling.classList.toggle('changed', !isDef);
			if (!isDef && !['file', 'seq', 'list'].includes(f.type)) changed[f.g] = (changed[f.g] || 0) + 1;
		}
		for (const b of el.querySelectorAll('[data-grp]')) { const n = changed[b.dataset.grp] || 0; b.textContent = n ? `${n} changed` : ''; b.classList.toggle('hidden', !n); }
	}
	return { sync };
}

const api = { GROUPS, build, setLine, toText };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuOptions = api;
})(typeof self !== 'undefined' ? self : globalThis);
