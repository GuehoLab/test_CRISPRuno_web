// SPDX-License-Identifier: LicenseRef-CRISPRuno — port of CRISPRunoCore.parse_settings (CRISPRuno 0.1.22, (c) MGH)
'use strict';
/* parse(text, { root, files }) -> { settings, warnings, used }
 *   text: a CRISPRuno settings file (key<TAB>value per line, '#' comments) — the command-line options of the
 *   original can be written as lines of the same file. File paths are reduced to their file names, which must be
 *   among the files selected on the page (FASTQ) — `genome` names the reference FASTA (optional here) and
 *   `aligner_genome` is ignored (the minibwa index is chosen on the page).
 *   Removed in the web port: Cas-OFFinder (PAM / casoffinder_num_mismatches are accepted and ignored),
 *   the *_command settings (accepted, ignored).
 */
(function (root) {
const C = root.CuCore || require('./core.js');

const DEFAULTS = [ // [key, default, type] in the order parse_settings fills settings{} (= order of <root>.settingsUsed.txt)
	['root', null, 'str'], ['debug', false, 'bool'], ['keep_intermediate', false, 'bool'], ['write_discarded_read_info', false, 'bool'], ['suppress_plots', false, 'bool'],
	['additional_cut_site_file', null, 'path'], ['cut_region_annotation_file', null, 'path'], ['cut_classification_annotations', [], 'annos'], ['PAM', null, 'str'],
	['guide_sequences', [], 'guides'], ['casoffinder_num_mismatches', 0, 'int'], ['cleavage_offset', -3, 'int'], ['primer_seq', null, 'str'], ['primer_in_r2', false, 'bool'],
	['min_primer_aln_score', 40, 'int'], ['allow_indels_in_origin_aln', false, 'bool'], ['min_primer_length', 30, 'int'], ['min_read_length', 30, 'int'],
	['transposase_adapter_seq', 'CTGTCTCTTATACACATCTGACGCTGCCGACGA', 'str'], ['min_alignment_quality_score', 30, 'int'], ['min_r2_alignment_quality_score', 30, 'int'],
	['arm_min_matched_start_bases', 10, 'int'], ['arm_max_clipped_bases', 0, 'int'], ['ignore_n', false, 'bool'], ['suppress_poor_alignment_filter', false, 'bool'],
	['guide_homology_max_gaps', 2, 'int'], ['guide_homology_max_mismatches', 5, 'int'], ['analyze_multimap_assignments', false, 'bool'], ['max_multimap_alignments', 5, 'int'],
	['run_crispresso_on_novel_sites', false, 'bool'], ['crispresso_min_count', 50, 'int'], ['crispresso_max_indel_size', 50, 'int'], ['crispresso_min_aln_score', 20, 'int'],
	['crispresso_quant_window_size', 1, 'int'], ['cutadapt_command', 'cutadapt', 'str'], ['samtools_command', 'samtools', 'str'], ['aligner_command', 'bowtie2', 'str'],
	['crispresso_command', 'CRISPResso', 'str'], ['casoffinder_command', 'cas-offinder', 'str'], ['n_processes', '1', 'nproc'],
	['novel_cut_merge_distance', 50, 'int'], ['known_cut_merge_distance', 50, 'int'], ['origin_cut_merge_distance', 10000, 'int'], ['suppress_homology_detection', false, 'bool'],
	['dedup_input_on_UMI', false, 'bool'], ['suppress_dedup_on_aln_pos_and_UMI_filter', false, 'bool'], ['dedup_by_final_cut_assignment_and_UMI', false, 'bool'],
	['umi_regex', null, 'str'], ['min_umi_seen_to_keep_read', 0, 'int'], ['write_UMI_counts', false, 'bool'], ['r1_r2_support_max_distance', 10000, 'int'],
	['suppress_r2_support_filter', false, 'bool'], ['fastq_r1', null, 'fastq'], ['fastq_r2', null, 'fastqopt'], ['fastq_umi', null, 'fastqopt'], ['genome', null, 'genome'], ['aligner_genome', null, 'str'],
];
// options of the command line that are not settings-file keys of the original but are accepted here as well
const CLI_ONLY = new Set(['short_indel_length_cutoff']);

function pyInt(key, v) {
	const s = String(v).trim();
	if (!/^[+-]?\d+(_\d+)*$/.test(s)) throw new Error(`ValueError: invalid literal for int() with base 10: '${v}' (setting ${key})`);
	return parseInt(s.replace(/_/g, ''), 10);
}
function splitList(val) { // guide_sequences / cut_classification_annotations as parse_settings splits them
	let arr = val.split(' ');
	if (val.includes('[')) arr = val.replace(/^[\[\] ]+|[\[\] ]+$/g, '').replace(/"/g, '').replace(/'/g, '').split(',').map((x) => C.strip(x));
	const out = [];
	for (const a of arr) { if (a.includes(',')) out.push(...a.split(',')); else if (a.includes(' ')) out.push(...a.split(' ')); else out.push(a); }
	return out.map((x) => x.replace(/^[\[\] '"]+|[\[\] '"]+$/g, ''));
}

function parse(text, opts = {}) {
	const files = opts.files || {}; // basename -> File
	const warnings = [];
	const args = new Map();
	for (const line of String(text).split('\n')) {
		if (line.startsWith('#')) continue;
		const els = line.split('#')[0].replace(/\r$/, '').split('\t');
		if (els.length === 1 && els[0] === '') continue;
		if (els.length < 2) throw new Error(`Cannot parse line "${line}"\nA tab must separate the key and value in the settings file`);
		const key = C.strip(els[0]), val = C.strip(els[1]);
		if (args.has(key) && args.get(key) !== val) throw new Error(`The parameter ${key} is provided twice`);
		args.set(key, val);
	}
	const S = {};
	for (const [key, def, type] of DEFAULTS) {
		const has = args.has(key), v = has ? args.get(key) : undefined;
		args.delete(key);
		if (type === 'bool') S[key] = has ? v.toLowerCase() === 'true' : def;
		else if (type === 'int') S[key] = has ? pyInt(key, v) : def;
		else if (type === 'str') S[key] = has ? v : def;
		else if (type === 'nproc') S[key] = has ? v : def;
		else if (type === 'path') S[key] = has ? (v === 'None' ? null : v) : null;
		else if (type === 'annos') {
			S[key] = [];
			if (has) for (let a of splitList(v)) {
				if (a === '') continue;
				if (a.split(':').length !== 4) throw new Error(`Error: cut classification annotations must contain 4 elements (e.g. chr1:234:left:my_annotation). Please check the value "${a}"`);
				S[key].push(a.replace(/__/g, ' '));
			}
		} else if (type === 'guides') {
			S[key] = [];
			if (has) for (const g of splitList(v)) {
				if (g === '') continue;
				const w = C.findWrongNt(g);
				if (w.length) throw new Error(`The guide sequence ${g} contains bad characters: "${w.join('", "')}"`);
				S[key].push(g);
			}
		} else if (type === 'fastq' || type === 'fastqopt') {
			let p = has ? v : null;
			if (type === 'fastqopt' && p && p.toLowerCase().includes('none')) p = null;
			S[key] = p;
		} else if (type === 'genome') S[key] = has ? v : null;
	}
	for (const k of CLI_ONLY) if (args.has(k)) { warnings.push(`${k} is a command-line option that CRISPRuno does not use (ignored)`); args.delete(k); }
	if (args.size) throw new Error('Unused keys in settings file: "' + [...args.keys()].join('", "') + '"');

	S.root = opts.root || (S.root ? C.base(S.root) : 'CRISPRuno');
	if (S.n_processes === 'max') S.n_processes = opts.cores || 1;
	else { const n = String(S.n_processes).trim(); if (!/^[+-]?\d+$/.test(n)) throw new Error(`Error: n_processes must be specified as and integer or "max" (current setting: ${S.n_processes})`); S.n_processes = parseInt(n, 10); }
	// files: only the names are used; they must be among the selected files
	const need = (key, p, required) => {
		if (!p) { if (required) throw new Error(`Error: ${key} file must be provided (--${key})`); return null; }
		const b = C.base(p);
		if (!files[b]) { if (opts.lenient) warnings.push(`${key}: ${b} is not among the selected files (needed by the steps that read it)`); else throw new Error(`Error: ${key} file ${b} is not among the selected files`); }
		return b;
	};
	S.fastq_r1 = need('fastq_r1', S.fastq_r1, true);
	S.fastq_r2 = need('fastq_r2', S.fastq_r2, false);
	S.fastq_umi = need('fastq_umi', S.fastq_umi, false);
	if (S.additional_cut_site_file) { need('additional_cut_site_file', S.additional_cut_site_file, true); S.additional_cut_site_file = C.base(S.additional_cut_site_file); warnings.push('additional_cut_site_file: CRISPRuno 0.1.22 reads this setting but does not pass it on to prep_input, so it has no effect (same here)'); }
	if (S.cut_region_annotation_file) S.cut_region_annotation_file = need('cut_region_annotation_file', S.cut_region_annotation_file, true);
	if (!S.primer_seq) throw new Error('Error: the primer sequences used must be provided (--primer_seq)');
	const w = C.findWrongNt(S.primer_seq);
	if (w.length) throw new Error(`The primer sequence ${S.primer_seq} contains bad characters: "${w.join('", "')}"`);
	if (S.PAM !== null || S.casoffinder_num_mismatches > 0) warnings.push('PAM / casoffinder_num_mismatches: the Cas-OFFinder search is not part of the web port and is skipped (only guide-aligned cut sites, homology and novel sites are used)');
	if (S.genome) S.genome = C.base(S.genome);
	return { settings: S, warnings };
}

// <root>.settingsUsed.txt (str() of each value, Python list repr for the lists)
function pyRepr(v) {
	if (Array.isArray(v)) return '[' + v.map((x) => "'" + String(x).replace(/\\/g, '\\\\').replace(/'/g, "\\'") + "'").join(', ') + ']';
	return C.str(v);
}
function settingsUsed(S) { return DEFAULTS.map(([k]) => `${k}\t${pyRepr(S[k])}\n`).join(''); }

const TEMPLATE = `# CRISPRuno settings (tab-separated key and value; file paths: only the file names are used)
fastq_r1\tsample_R1.fastq.gz
fastq_r2\tsample_R2.fastq.gz
fastq_umi\tsample_UMI.fastq.gz
genome\thg38.fa
primer_seq\tAGGGCACACCCTGTGATCTTGTGGG
guide_sequences\tTGCTTGGTCGGCACTGATAG
min_primer_length\t25
`;

const api = { parse, settingsUsed, DEFAULTS, TEMPLATE };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuSettings = api;
})(typeof self !== 'undefined' ? self : globalThis);
