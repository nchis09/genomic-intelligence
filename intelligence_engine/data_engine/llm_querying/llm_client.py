"""
llm_client.py — Local Ollama LLM client for epidemiological context queries.

Only Ollama is supported for data privacy. No data leaves the machine.
Default model: Qwen 2.5 14B via Ollama.

Usage:
    from intelligence_engine.data_engine.llm_querying.llm_client import LLMClient

    client = LLMClient()  # auto-detects local Ollama
    response = client.query_epidemiology(prompt)
    # response is a parsed dict (JSON mode enabled)

    report = client.write_report(intelligence_object_str)
    # report is a string (narrative text)
"""

import json
import logging
import os
from typing import Optional

log = logging.getLogger(__name__)

# ─── Provider configs ─────────────────────────────────────────────────────────

OLLAMA_URL = os.environ.get("OLLAMA_URL", "http://localhost:11434")
OLLAMA_MODEL = os.environ.get("OLLAMA_MODEL", "qwen2.5:7b")
# How long Ollama keeps the model in memory after the last request.
# Default 30m avoids a full reload on every short engine run. Set to 0 or "5m" as needed.
OLLAMA_KEEP_ALIVE = os.environ.get("OLLAMA_KEEP_ALIVE", "30m")

class LLMClient:
    """Local Ollama LLM client.

    Connects to a local Ollama server for privacy. If Ollama is not available,
    callers should fall back to deterministic mode or raise an error.
    """

    def __init__(self) -> None:
        self.provider = self._resolve_provider()
        log.info(f"LLM client initialised: provider={self.provider}, "
                 f"model={OLLAMA_MODEL}")

    def _resolve_provider(self) -> str:
        """Confirm that Ollama is available."""
        if self._ollama_available():
            return "ollama"
        raise RuntimeError(
            "Ollama is not available. Install Ollama, run `ollama pull qwen2.5:14b`, "
            "and start `ollama serve` before enabling the LLM."
        )

    def _ollama_available(self) -> bool:
        """Check if Ollama is running and the model is available."""
        import urllib.request
        use_llm_env = os.environ.get("PGIRL_USE_LLM", "true").lower()
        if use_llm_env in ("false", "0", "no", "off"):
            log.info("PGIRL_USE_LLM=%s; skipping Ollama availability check.", use_llm_env)
            return False
        try:
            req = urllib.request.Request(f"{OLLAMA_URL}/api/tags", method="GET")
            with urllib.request.urlopen(req, timeout=3) as resp:
                data = json.loads(resp.read())
                models = [m.get("name", "") for m in data.get("models", [])]
                # Check if our model (or a prefix match) is available
                for m in models:
                    if OLLAMA_MODEL in m or m in OLLAMA_MODEL:
                        return True
                log.warning(f"Ollama running but model '{OLLAMA_MODEL}' not found. "
                            f"Available: {models}")
                return False
        except Exception:
            return False

    def _get_model(self) -> str:
        return OLLAMA_MODEL

    def _get_url(self) -> str:
        return f"{OLLAMA_URL}/api/chat"

    def _get_headers(self) -> dict:
        return {"Content-Type": "application/json"}

    def _build_payload(self, messages: list, temperature: float,
                       max_tokens: int, json_mode: bool) -> dict:
        """Build request payload for the active provider."""
        payload = {
            "model": self._get_model(),
            "messages": messages,
            "temperature": temperature,
            "max_tokens": max_tokens,
        }
        if self.provider == "ollama":
            payload["keep_alive"] = OLLAMA_KEEP_ALIVE
        if json_mode:
            if self.provider == "ollama":
                payload["format"] = "json"
            else:
                payload["response_format"] = {"type": "json_object"}
        return payload

    def _call_api(self, payload: dict, timeout: Optional[int] = None) -> str:
        """Make the HTTP request and return the raw text response."""
        import urllib.request
        # Disable streaming for Ollama — we want a single JSON response
        if self.provider == "ollama":
            payload["stream"] = False
        body = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            self._get_url(),
            data=body,
            headers=self._get_headers(),
            method="POST",
        )
        if timeout is None:
            timeout = 300
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read()

        data = json.loads(raw)
        return data.get("message", {}).get("content", "")

    def _call(self, system_prompt: str, user_prompt: str,
              temperature: float = 0.3, max_tokens: int = 4096,
              json_mode: bool = False, timeout: Optional[int] = None) -> str:
        """Call the LLM and return raw text response.

        Args:
            system_prompt: System instruction (role, constraints)
            user_prompt: The actual questions/data
            temperature: 0.3 for factual queries, 0.7 for report writing
            max_tokens: Response token limit
            json_mode: If True, force JSON output
            timeout: Optional override for the HTTP request timeout (seconds).
                Useful for long-form generations (e.g. full reports) that
                exceed the default provider timeout on local hardware.
        """
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ]
        payload = self._build_payload(messages, temperature, max_tokens, json_mode)

        return self._call_api(payload, timeout=timeout)

    def query_epidemiology(self, user_prompt: str) -> dict:
        """Call LLM for epidemiological context. Returns parsed JSON dict.

        The LLM is instructed to return structured JSON with numbers where
        possible and null for unknown data. Temperature is low (0.3) to
        minimise hallucination.

        Args:
            user_prompt: The assembled prompt with bioinformatics summary,
                         DB results, and conditional questions.

        Returns:
            Parsed JSON dict with epidemiological context.
        """
        system_prompt = (
            "You are an infectious disease epidemiologist acting as a strict data extraction assistant. "
            "Use ONLY the fetched epidemiological text provided by the user. "
            "Do NOT use your training data, memory, prior knowledge, or general world facts to answer. "
            "Before using any piece of information, verify that the source is an official or trusted public-health "
            "page (e.g. WHO, CDC, Africa CDC, ECDC, GOARN, national public-health agencies, or recognised international organisations). "
            "Ignore or discard claims from unofficial, unknown, or low-credibility sources. "
            "Answer only the epidemiological questions that are explicitly listed in the prompt. "
            "The fetched data is ordered by source credibility. Prefer WHO/CDC/official agency reports first. "
            "Avoid repeating the same outbreak or event; pick the most important and most complete entries for each question. "
            "If multiple sources report the same event, choose the one with the most complete numbers and cite its source URL. "
            "If the provided text does not contain enough information for a question, return an empty array [] for that key. "
            "Return STRUCTURED JSON with numbers where possible. "
            "Cite source URLs from the fetched data where possible. "
            "Do NOT make up data, infer beyond the provided text, or fill gaps with prior knowledge. "
            "Return ONLY valid JSON, no prose before or after."
        )

        raw = self._call(
            system_prompt=system_prompt,
            user_prompt=user_prompt,
            temperature=0.3,
            max_tokens=4096,
            json_mode=True,
        )

        return self._parse_json_response(raw)

    def write_report(self, intelligence_object_str: str) -> str:
        """Call LLM to write a narrative report from structured intelligence.

        This is the second LLM call — it only summarises the structured
        findings, does not do analysis. Temperature is higher (0.7) for
        more natural prose.

        Args:
            intelligence_object_str: JSON string of the combined intelligence object.

        Returns:
            Narrative report text (GOARN field report format).
        """
        system_prompt = (
            "You are a scientific writer for a public health agency. "
            "Given the structured intelligence object below, write a genomic intelligence brief "
            "in the format of a GOARN field report. "
            "Use the biological facts, epidemiological context, risk assessment, and limitations. "
            "Do NOT include a recommended actions or response recommendations section. "
            "Include section headers (S1-S8): S1 Sample Information, S2 Genomic Characteristics, "
            "S3 Variant Summary, S4 Historical Lineage Context, S5 Epidemiological Overview, "
            "S6 Risk Assessment, S7 Limitations, S8 References. "
            "Do not add information not present in the data. "
            "Be precise and factual."
        )

        return self._call(
            system_prompt=system_prompt,
            user_prompt=intelligence_object_str,
            temperature=0.7,
            max_tokens=4096,
            json_mode=False,
        )

    @staticmethod
    def _grounding_and_reasoning_block() -> str:
        """Shared grounding/reasoning instructions used by both the Brief and
        the full Report synthesis prompts. Kept in one place so the two
        document types cannot silently drift apart in their evidence rules."""
        return (
            "WHAT REASONING MEANS HERE: You SHOULD use your full reasoning and writing ability to connect "
            "facts that are logically related, notice patterns across evidence streams, and produce a novel, "
            "well-organized synthesis — this is expected and valuable. What you must NOT do is supply any "
            "pathogen-specific FACT (a number, name, date, mutation, outbreak, lineage, location, citation) "
            "that is not literally present in the evidence JSON below. Reasoning/logic is yours to apply; "
            "facts must always come from the JSON, never from memory of this or any other pathogen.\n"
            "\n\n"
            "STRICT GROUNDING RULES:\n"
            "1. Use ONLY the provided evidence base for facts. Do NOT use your training data, memory, prior "
            "knowledge, or general world facts about this pathogen/lineage/species to add information not "
            "present in the evidence, even if you believe it to be true. If the evidence does not cover "
            "something, say so explicitly as a knowledge gap rather than filling it from memory.\n"
            "2. Do NOT perform new statistical analyses, invent numbers, or change units/values. "
            "Quote the exact numbers from the evidence. If a value is approximate, use the wording "
            "provided in the evidence.\n"
            "3. Do NOT assign a risk tier (e.g. MONITOR / INVESTIGATE / HIGH PRIORITY), predict or "
            "recommend public-health actions, response measures, or outbreak investigation steps. Do NOT "
            "tell the reader what they 'should', 'must', or 'need to' do. This is an evidence product, not "
            "an operational plan or a risk classification — the epidemiologist makes the risk assessment; "
            "you provide the specific, connected, and convergence-tested facts they need to make it.\n"
            "4. Be maximally specific, never ambiguous. Do not write vague statements like 'there were "
            "many outbreaks' or 'the region has had past cases'. Instead state exactly what the evidence "
            "says: how many, where (country/admin region), when (dates/years), who (host/species/cases/"
            "deaths), and cite it. Every factual claim should answer as many of what/where/when/who/how-many "
            "as the evidence provides.\n"
            "\n\n"
            "EVIDENCE CONVERGENCE — this is central to the reasoning style: for every major conclusion, "
            "explicitly state which INDEPENDENT evidence streams (genomic/bioinformatics, phylogenetic, "
            "molecular-epidemiological, historical-outbreak, cross-evidence statistics) agree with each "
            "other, and which conflict or are silent. A conclusion supported by 3 independent streams should "
            "read differently from one resting on a single source — say so explicitly (e.g. 'supported "
            "independently by lineage assignment, phylogenetic placement, and genetic relatedness "
            "analysis...' vs 'based solely on ...'). Do not present findings as a flat sequential list; "
            "show the reasoning chain that connects them.\n"
            "\n\n"
            "QUANTITATIVE EVIDENCE RULES:\n"
            "- Every substantive statement MUST be supported by explicit quantitative evidence from "
            "  the JSON: genome completeness, depth, SNP distances, percent identity, lineage prevalence, "
            "  mutation frequency/counts, temporal persistence values, molecular-clock estimates, confidence "
            "  scores, p-values, sample sizes, or literature references.\n"
            "- Do NOT use vague qualitative descriptors ('high', 'low', 'significant', 'notable', "
            "  'concerning', 'substantial', 'important', 'major', 'widespread', 'rare', 'several', 'various', "
            "  'multiple', 'primarily', 'limited', 'extensive', 'strong', 'weak') as a stand-in for a number. "
            "  State the number itself as the sentence, e.g. write 'the sequence reached 99.56% completeness "
            "  at a mean depth of 142x, with all genes covered except the 3' UTR' -- NOT 'the sequence was "
            "  high quality (99.56% completeness, 142x depth)'. The adjective is redundant once the number "
            "  is stated; drop the adjective entirely and let the number carry the sentence.\n"
            "- Quote exact values directly in the sentence, not parked in a parenthetical after a "
            "  qualitative claim -- the number IS the claim.\n"
            "\n\n"
            "NO INLINE CITATION TAGS:\n"
            "- Do NOT append inline citation tags such as '(evidence: xxx)', '(source: xxx)', or "
            "  '[PMID:xxx]' after claims. Write the finding directly as a stated fact in plain prose -- "
            "  the reader should never see an internal evidence-stream, table, or finding_type name.\n"
            "- The one exception is literature: when a specific publication supports a claim, cite it "
            "  naturally in-line as '(PMID: <id>)' immediately after that claim, the way a scientific "
            "  report would cite a reference -- not as a generic evidence tag.\n"
            "- Every fact must still be grounded in the evidence JSON below (see STRICT GROUNDING RULES); "
            "  citation tags are simply not shown to the reader. Full traceability is preserved separately "
            "  in context_used.json.\n"
            "- Do not invent facts to avoid a citation -- ground everything, just state it as prose.\n"
            "\n\n"
            "ACCESSION / IDENTIFIER RULE:\n"
            "- When mentioning any accession number, genome ID, or strain identifier, ALWAYS precede "
            "  it with the organism/pathogen/strain name it refers to. Never write a bare accession "
            "  like 'KC242791.1' without context. Instead write: 'Makona-GH-MAL-2014 (accession "
            "  KC242791.1)' so the reader knows what it is without cross-referencing.\n"
            "\n\n"
            "TRADITIONAL EPIDEMIOLOGY RULE — CRITICAL:\n"
            "- The evidence JSON contains outbreak-level epidemiological data (cases, deaths, CFR, "
            "  country, year) in multiple locations: epidemiological_query_summary.outbreaks, "
            "  epidemiological_summary.preview_blocks, genomic_analyses.comparative_outbreak_analysis, "
            "  and cross_evidence_statistics.\n"
            "- You MUST incorporate these traditional epidemiological metrics (case counts, death "
            "  counts, CFR percentages, outbreak years, affected countries) into the Historical and "
            "  Epidemiological Context sections. A public health officer reading the brief/report needs "
            "  to see the scale of impact of this pathogen/lineage, not just molecular data.\n"
            "- For each historical outbreak mentioned, state: country, year(s), cases, deaths, and CFR "
            "  if available — all from the evidence.\n"
            "- Distinguish molecular epidemiology (lineages, mutations, phylogeography) from traditional "
            "  epidemiology (outbreaks, cases, deaths, CFR, surveillance coverage). Both must be "
            "  present in the epidemiological sections.\n"
            "\n"
            "LINEAGE-SPECIFIC OUTBREAK LINKING — CRITICAL:\n"
            "- When reporting historical outbreaks, you MUST distinguish outbreaks caused by THIS "
            "  specific lineage/species from outbreaks caused by OTHER lineages or species of the "
            "  same pathogen family. An epidemiologist needs to know the historical impact of THIS "
            "  lineage to assess potential damage from the current detection.\n"
            "- Link the molecular identity (lineage, clade, species) to the epidemiological data: "
            "  'This lineage (EBOV-Ebov-2013) caused X cases and Y deaths in Z countries during "
            "  <years>' rather than generically listing any Ebola outbreak in the region.\n"
            "- If outbreaks in the same country were caused by a different lineage or species, state "
            "  this explicitly: 'The 2025 Uganda outbreak was caused by a different species/lineage, "
            "  not EBOV-Ebov-2013' so the reader does not conflate them.\n"
            "- The comparative_outbreak_analysis in the evidence ranks outbreaks by similarity to the "
            "  current sample — use this to identify which outbreaks are most relevant to THIS "
            "  lineage, not just geographically coincidental.\n"
            "\n"
            "VARIANT-TO-OUTBREAK-IMPACT LINKING — CRITICAL:\n"
            "- If the evidence JSON contains 'genomic_links' entries connecting a detected mutation or "
            "  lineage to a specific outbreak_id, use them to state exactly which outbreak(s) carried "
            "  this mutation, then pull that outbreak's cases, deaths, CFR, country, and year from "
            "  'outbreaks' or 'molecular_epidemiology' and report them together as one connected fact.\n"
            "- If the evidence JSON contains a 'demographics' section (age_group, sex, occupation, "
            "  population_affected, risk_group, exposure_history, setting, case_count), state who was "
            "  affected — age group, sex, occupation/risk group, exposure setting — for the outbreak(s) "
            "  relevant to this lineage/variant. If the evidence does not tie a demographic breakdown to "
            "  a specific outbreak, state the demographic finding and the outbreak/country/year it came "
            "  from side by side so the reader can judge the link themselves.\n"
            "- If no demographic breakdown (age/sex/case counts by group) exists in the evidence for "
            "  this lineage's outbreaks, state explicitly that this breakdown is not available in the "
            "  evidence — do not omit the question or answer it from memory.\n"
            "\n"
            "SURVEILLANCE INTERPRETATION RULE:\n"
            "- Do NOT assume that an absence of genome sequences means a surveillance gap or failure. "
            "  If there were no reported cases of a pathogen in a region during a period, the absence "
            "  of sequences is expected and normal, not evidence of poor monitoring.\n"
            "- Only describe a surveillance gap when the evidence explicitly indicates missing "
            "  coverage despite known ongoing transmission. If cases were not reported, state that "
            "  no cases were reported, rather than attributing the absence to monitoring failures.\n"
        )

    @staticmethod
    def _figures_instruction_block(figures: Optional[list]) -> str:
        """Build the figure-reference instructions from the exact figures the
        pipeline generated, so the model cannot invent figure numbers/names."""
        if not figures:
            return (
                "FIGURES: No figures were generated for this sample. Omit the VISUALISATIONS section "
                "entirely (do not write a VISUALISATIONS heading with no content).\n"
            )
        lines = [
            "FIGURES: The following figures were actually generated for this sample. Use ONLY these, "
            "with these exact numbers and titles — do not invent, renumber, or add any other figure:",
        ]
        for f in figures:
            lines.append(f"  Figure {f.get('figure_number')} — {f.get('title')} (file: {f.get('filename')})")
        lines.append(
            "In the VISUALISATIONS section, reference each figure exactly as: "
            "'Ref: Figure <N> – <Title>' on its own line, followed by the filename and a one-line "
            "description of what it shows grounded in the evidence (not invented commentary).\n"
        )
        return "\n".join(lines)

    def synthesize_intelligence_brief(
        self,
        evidence_context_str: str,
        brief_id: str,
        report_id: str,
        sample_id: str,
        generated_at: str,
        figures: Optional[list] = None,
    ) -> str:
        """Generate the Genomic Intelligence Brief as a plain-text ASCII file.

        Combines structured header fields (PATHOGEN IDENTIFICATION, SAMPLE
        CONTEXT) with narrative reasoning sections (Executive Assessment through
        Overall Genomic Intelligence Assessment), using '=' and '-' dividers
        for visual clarity. No Markdown syntax.
        """
        system_prompt = (
            "You are a senior genomic intelligence analyst drafting a GENOMIC INTELLIGENCE BRIEF — a "
            "concise executive snapshot for a public health specialist / epidemiologist who will use it "
            "as the factual basis for their own rapid risk assessment. It must be readable in under five "
            "minutes but still show the reasoning connections between evidence streams, in the spirit of "
            "a WHO Rapid Risk Assessment executive summary.\n\n"
            + self._grounding_and_reasoning_block()
            + "\n\n"
            "REASONING STYLE:\n"
            "- The EXECUTIVE ASSESSMENT is the single most important section: it must read as connected "
            "  reasoning that links genomic findings (mutations/phenotypes), molecular epidemiology "
            "  (lineage history), classical epidemiology (outbreaks/cases/deaths/demographics), and "
            "  intervention availability (vaccines/therapeutics) into one coherent narrative -- this is "
            "  the genomic intelligence linkage capability the reader is paying for, not a restatement "
            "  of module outputs in isolation.\n"
            "- Introduce quantitative values only when they directly support an interpretation. Explain "
            "  what the number means, not just what it is.\n"
            "- Integrate across domains. Show how lineage assignment, phylogenetic placement, genetic "
            "  relatedness, molecular clock, and mutation profile converge or diverge. Identify which "
            "  conclusions are supported by multiple independent evidence streams and which rest on a "
            "  single source.\n"
            "- Distinguish types of uncertainty: (a) scientific uncertainty (uncharacterised mutations, "
            "  unknown phenotypes), (b) surveillance uncertainty (missing recent genomes, missing Rt/R0 "
            "  estimates), and (c) analytical limitation (insufficient co-occurrence data, database gaps).\n"
            "\n"
            "QUANTIFICATION RULE — CRITICAL:\n"
            "- State the exact quantitative value as the sentence itself, not parked in a bracket after a "
            "qualitative adjective. The reader has no time to cross-check with the full report.\n"
            "- Examples of correct quantification (structure only — use the ACTUAL values from the "
            "  evidence JSON, never these placeholder examples):\n"
            "  'the sequence reached <completeness>% completeness at <depth>x mean depth'\n"
            "  '<identity>% identical to <reference strain>, differing by <N> SNPs'\n"
            "  'placed with <bootstrap>% bootstrap support'\n"
            "  'the lineage comprises <count> curated genomes spanning <year range>'\n"
            "  'the outbreak recorded <cases> cases and <deaths> deaths (<CFR>% CFR)'\n"
            "- Do not write qualitative adjectives like 'high', 'close', 'similar', 'strong', "
            "'well-supported', 'large', or 'significant' at all; replace them with the number itself, "
            "e.g. write '98% bootstrap support' rather than 'well-supported (98% bootstrap)'.\n"
            "- This applies to ALL sections: PATHOGEN IDENTIFICATION, SAMPLE CONTEXT, EXECUTIVE "
            "ASSESSMENT, KEY GENOMIC FINDINGS, and EPIDEMIOLOGICAL CONTEXT.\n"
            "\n\n"
            + self._figures_instruction_block(figures)
            + "\n\n"
            "OUTPUT FORMAT — you MUST reproduce this exact plain-text template, replacing the placeholder "
            "text in <angle brackets> with content grounded in the evidence JSON. Keep the '=' and '-' "
            "divider lines EXACTLY as shown (80 characters). Do not use Markdown syntax (no #, no **, no "
            "Markdown bullets); use plain '-' or indented dashes as in the template. Do not add any section "
            "not listed below. If a field has no supporting evidence, write 'Not available in evidence' "
            "rather than omitting the line or inventing a value.\n\n"
            "================================================================================\n"
            "                         GENOMIC INTELLIGENCE BRIEF\n"
            "================================================================================\n"
            f"Brief ID:      {brief_id}\n"
            f"Report ID:     {report_id}\n"
            f"Generated:     {generated_at}\n"
            f"Sample ID:     {sample_id}\n"
            "================================================================================\n\n"
            "PATHOGEN IDENTIFICATION\n"
            "--------------------------------------------------------------------------------\n"
            "Species:          <species>\n"
            "Lineage / Clade:  <lineage> | <clade>\n"
            "Confidence:       <confidence level, grounded in evidence>\n"
            "Closest genome:   <SNP distance, identity %, and accession/strain from the evidence>\n"
            "Lineage in DB:    <total curated genomes, first-last detection dates, countries reported>\n\n"
            "SAMPLE CONTEXT\n"
            "--------------------------------------------------------------------------------\n"
            "Location:         <country/admin1/admin2 from evidence>\n"
            "Collection date:  <date>\n"
            "Source / host:    <sampling source, host>\n"
            "Genome quality:   <completeness %, depth, quality flag, missing regions>\n\n"
            "EXECUTIVE ASSESSMENT\n"
            "--------------------------------------------------------------------------------\n"
            "  <A dense, multi-paragraph synthesis (roughly half to three-quarters of an A4 page) "
            "written as connected reasoning, not a checklist -- this is where the genomic intelligence "
            "linkage across ALL evidence sources (bioinformatics output, database query output, "
            "evidence-integration analyses, and any user-provided sample context) must be visible. "
            "It must, in prose grounded strictly in the evidence JSON below:\n"
            "  - State what was found: species/lineage/clade, confidence, closest reference and closest "
            "outbreak/contextual genome with their exact identity%/SNP distances, and which independent "
            "evidence streams (phylogenetic placement, genetic relatedness, lineage assignment) converge "
            "to support the identification.\n"
            "  - Name ONLY the mutations that carry curated phenotype or literature evidence in the "
            "evidence JSON (state how many total variants were detected and how many of those have no "
            "curated evidence, without listing the uncharacterised ones by name here). For each mutation "
            "with evidence, state its exact effect category, evidence strength, frequency in the curated "
            "database (genome_count/total_genomes), and the specific biological/clinical outcome that "
            "phenotype implies (e.g. immune escape, altered vaccine effectiveness, increased "
            "transmissibility), exactly as stated in the evidence.\n"
            "  - State the molecular epidemiology finding: where and when this lineage was last detected, "
            "in which countries, and how many curated genomes represent it.\n"
            "  - State the classical epidemiology finding for outbreaks caused by this species/lineage: "
            "cases, deaths, CFR, countries, and years (from 'outbreaks', 'molecular_epidemiology', "
            "'comparative_outbreak_analysis'), AND, where available in the evidence, who was affected -- "
            "age group, sex, occupation/risk group, exposure setting -- from 'demographics', linked via "
            "'genomic_links' to this lineage where the evidence supports it. If no demographic breakdown "
            "exists for this lineage's outbreaks, state that explicitly.\n"
            "  - State vaccine/therapeutic/intervention availability and effectiveness exactly as recorded "
            "in 'vaccines', 'interventions', and 'therapeutics'.\n"
            "  - State the phylogeographic/introduction assessment (inferred origin, its probability/"
            "support values, dissemination routes, temporal gap in days) with the exact numbers.\n"
            "  - Close with one sentence stating what remains uncertain and why.\n"
            "  Do not state a conclusion, risk tier, or recommendation anywhere in this section. Every "
            "sentence must carry a specific number, date, place, or named source from the evidence -- no "
            "unsupported adjectives, no restating a fact as merely 'high' or 'low' without its value.>\n\n"
            "KEY GENOMIC FINDINGS\n"
            "--------------------------------------------------------------------------------\n"
            "  <one block per mutation that HAS curated evidence in the JSON (phenotype associations, "
            "hotspot status, or literature): mutation (gene:change), domain, known/novel status, effect "
            "category, evidence strength, frequency (genome_count/total_genomes), first/last seen dates, "
            "countries seen, and the specific literature finding(s) stated as fact with '(PMID: <id>)' "
            "where given. For mutations with no curated evidence, list them by name only, one line: "
            "'No curated evidence available for: <mutation list>'.>\n\n"
            "EPIDEMIOLOGICAL CONTEXT\n"
            "--------------------------------------------------------------------------------\n"
            "Historical cases/deaths/CFR:    <exact numbers for this species/lineage's outbreaks>\n"
            "Affected populations:           <age/sex/occupation/exposure from 'demographics', or "
            "'Not available in evidence'>\n"
            "Licensed vaccine(s):            <product, status, effectiveness from 'vaccines', or "
            "'Not available in evidence'>\n"
            "Other interventions:            <from 'interventions'/'therapeutics', or 'Not available "
            "in evidence'>\n"
            "Transmission route:             <from 'transmission' evidence>\n"
            "Introduction assessment:        <new introduction / re-emergence / local persistence, with "
            "the exact scores from 'introduction_scenarios' and the inferred origin>\n\n"
            "CRITICAL KNOWLEDGE GAPS\n"
            "--------------------------------------------------------------------------------\n"
            "  <3-6 explicit, categorised gaps (scientific / surveillance / analytical), each stating "
            "exactly what data is missing, from 'knowledge_gaps' and any other explicit gaps visible in "
            "the evidence>\n\n"
            "REFERENCES\n"
            "--------------------------------------------------------------------------------\n"
            "  <numbered list of literature/database references actually used above, from 'references' "
            "in the evidence>\n\n"
            + (
                "VISUALISATIONS\n"
                "--------------------------------------------------------------------------------\n"
                "  <'Ref: Figure <N> – <Title>' lines per the FIGURES instructions above, each followed "
                "by the filename and a one-line grounded description; omit section if none>\n\n"
                if figures else ""
            )
            + "================================================================================\n"
            "  This brief is auto-generated. It MUST be reviewed and approved by a\n"
            "  qualified public health genomicist before dissemination or action.\n"
            f"  Full evidence chain: see Genomic Intelligence Report {report_id}\n"
            "================================================================================\n"
        )

        hard_constraints_reminder = self._hard_constraints_reminder(ascii_mode=True)

        raw = self._call(
            system_prompt=system_prompt,
            user_prompt=evidence_context_str + self._user_preface() + hard_constraints_reminder,
            temperature=0.2,
            max_tokens=7168,
            json_mode=False,
            timeout=480,
        )
        return self._enforce_ascii_contract(raw)

    def synthesize_intelligence_report(
        self,
        evidence_context_str: str,
        report_id: str,
        brief_id: str,
        sample_id: str,
        generated_at: str,
        figures: Optional[list] = None,
    ) -> str:
        """Generate Part 2 — the full Genomic Intelligence Report — following
        the project's fixed-width ASCII template (see README.md, "Part 2:
        Genomic Intelligence Report"), with a dedicated Evidence Convergence
        section and without any risk-tier classification subsection.
        """
        system_prompt = (
            "You are a senior genomic intelligence analyst drafting the full GENOMIC INTELLIGENCE REPORT — "
            "the comprehensive evidence document behind a Genomic Intelligence Brief, for analysts/experts "
            "who need the complete evidence chain to conduct their own risk assessment. Every section must "
            "read like connected reasoning (in the spirit of a WHO Rapid Risk Assessment), not a sequential "
            "dump of module outputs: explicitly show how independent evidence streams converge or conflict "
            "to support each conclusion.\n\n"
            + self._grounding_and_reasoning_block()
            + "\n\n"
            + self._figures_instruction_block(figures)
            + "\n\n"
            "QUANTIFICATION RULE: State exact quantitative values directly as the sentence itself, not "
            "parked in a bracket after a qualitative adjective, e.g. 'the sequence reached "
            "<completeness>% completeness at <depth>x mean depth' rather than 'high-quality sequence "
            "(<completeness>% completeness, <depth>x depth)'. Use the ACTUAL values from the evidence "
            "JSON, never placeholder or example values.\n\n"
            "OUTPUT FORMAT — reproduce this exact plain-text template, replacing placeholder text in "
            "<angle brackets> with content grounded in the evidence JSON. Keep the '=' and '-' divider "
            "lines EXACTLY as shown. Do NOT use Markdown syntax — no **, no #, no backticks, no "
            "Markdown bullets. Use 'Not available in evidence' for uncovered "
            "fields instead of omitting or inventing. Do not add a Risk Classification / escalation-trigger "
            "subsection, nor an 'Evidence-Based Considerations for Decision-Making' / recommendations-style "
            "section, anywhere in the report -- this is an evidence-synthesis document, not an operational "
            "one.\n\n"
            "================================================================================\n"
            "                      GENOMIC INTELLIGENCE REPORT\n"
            "================================================================================\n"
            f"Report ID:     {report_id}\n"
            f"Brief ID:      {brief_id}\n"
            f"Generated:     {generated_at}\n"
            f"Sample ID:     {sample_id}\n"
            "Prepared by:   Genomic Epidemic Intelligence System (automated)\n"
            "Review status: PENDING_EXPERT_REVIEW\n"
            "================================================================================\n\n"
            "SECTION 1: SAMPLE & SURVEILLANCE CONTEXT\n"
            "--------------------------------------------------------------------------------\n"
            "Sample ID:              <sample id>\n"
            "Sampling source:        <source, e.g. clinical (blood, RT-PCR positive)>\n"
            "Collection date:        <date>\n"
            "Country:                <country>\n"
            "Admin Level 1:          <admin1, or 'Not available in evidence'>\n"
            "Admin Level 2:          <admin2, or 'Not available in evidence'>\n"
            "Host:                   <host species/status>\n"
            "Sequencing platform:    <platform, or 'Not available in evidence'>\n"
            "Sequencing protocol:    <protocol, or 'Not available in evidence'>\n"
            "Submitting lab:         <lab, or 'Not available in evidence'>\n"
            "Genome completeness:    <completeness % (bp/expected bp)>\n"
            "Mean coverage depth:    <depth>x\n"
            "Missing regions:        <region, length, note; or 'None reported'>\n"
            "Genome quality flag:    <flag>\n\n"
            "SECTION 2: GENOMIC CHARACTERIZATION\n"
            "--------------------------------------------------------------------------------\n"
            "Family:                  <family>\n"
            "Genus:                   <genus>\n"
            "Species:                 <species>\n"
            "Lineage:                 <lineage>\n"
            "Clade:                   <clade>\n"
            "Closest reference:       <accession> (<name>)\n"
            "  - Identity:            <identity %>\n"
            "  - SNPs from reference: <count>\n"
            "Closest outbreak/contextual genome: <accession> (<name>)\n"
            "  - SNPs from closest:   <count>\n"
            "Phylogenetic placement:  <cluster/clade placement>\n"
            "  - Support:             <bootstrap/UFboot value>\n"
            "Genome length:           <length bp>\n"
            "Genes covered:           <list> (<n>/<total>)\n"
            "Missing regions:         <region, length bp>\n"
            "Confidence:              <HIGH/MODERATE/LOW per dimension: species, lineage, clade>\n\n"
            "SECTION 3: FUNCTIONAL GENOMICS — MUTATIONS & GENOMIC FEATURES\n"
            "--------------------------------------------------------------------------------\n"
            "Total amino acid substitutions: <count>\n\n"
            "  3a. MUTATIONS WITH CURATED EVIDENCE\n"
            "  ...................................\n"
            "  <one block per mutation that HAS curated phenotype/literature evidence, formatted as:\n"
            "    Mutation: <gene:change>\n"
            "    Type:     <amino acid substitution (gene, domain)>\n"
            "    Effect:   <effect category exactly as stated in evidence>\n"
            "    Evidence:\n"
            "      - <finding stated as fact> (PMID: <id>) [<evidence_strength>]\n"
            "      - <repeat for each literature/phenotype finding>\n"
            "    Confidence: <evidence_strength/confidence_score from evidence>\n"
            "    Frequency: <genome_count>/<total_genomes> (<pct>%) across <countries_seen>, "
            "first seen <date>, last seen <date>\n"
            "    Public health relevance:\n"
            "      - Transmissibility: <stated fact or 'Not available in evidence'>\n"
            "      - Virulence:        <stated fact or 'Not available in evidence'>\n"
            "      - Diagnostics:      <stated fact or 'Not available in evidence'>\n"
            "      - Therapeutics:     <stated fact or 'Not available in evidence'>\n"
            "      - Vaccine:          <stated fact or 'Not available in evidence'>\n"
            "    Outbreak linkage: <where 'genomic_links' ties this mutation/lineage to a specific "
            "outbreak_id, name the outbreak(s) with cases/deaths/CFR; otherwise 'Not available in "
            "evidence'>>\n\n"
            "  3b. NOVEL / UNCHARACTERISED MUTATIONS\n"
            "  ......................................\n"
            "  <list mutations with NO curated evidence: gene:change, domain, nucleotide change; note "
            "'No peer-reviewed evidence available' for each>\n\n"
            "  3c. MUTATIONS NOT DETECTED (NEGATIVE FINDINGS)\n"
            "  ..............................................\n"
            "  <any explicit negative findings in the evidence, e.g. known concerning mutations checked "
            "for but absent; otherwise 'None stated in evidence'>\n\n"
            "SECTION 4: MOLECULAR EPIDEMIOLOGY & OUTBREAK HISTORY\n"
            "--------------------------------------------------------------------------------\n"
            "  4a. HISTORICAL OCCURRENCE\n"
            "  ..........................\n"
            "  Species/lineage first documented: <year, location>\n"
            "  Total confirmed outbreaks:         <count>\n"
            "  Cumulative cases / deaths / CFR:   <numbers, from 'epidemic_dynamics'/'outbreaks'>\n\n"
            "  4b. LINEAGE DISTRIBUTION\n"
            "  .........................\n"
            "  Lineage <id>:\n"
            "    - First detected:  <date, location>\n"
            "    - Last detected:   <date, location>\n"
            "    - Countries:       <list>\n"
            "    - Total genomes:   <count>\n"
            "    - Associated outbreaks: <one line per outbreak: name/date/country, cases, CFR>\n\n"
            "  4c. GEOGRAPHIC CONTEXT\n"
            "  .......................\n"
            "  Endemic in:  <countries>\n"
            "  This sample's location: previous detections of this species/lineage here (or 'first "
            "detection'), proximity to nearest prior outbreak in km/days if stated.\n\n"
            "  4d. RESERVOIR & HOST RANGE\n"
            "  ............................\n"
            "  Known reservoirs:   <from evidence, or 'Not available in evidence'>\n"
            "  Accidental hosts:   <from evidence, or 'Not available in evidence'>\n"
            "  Transmission route: <from 'transmission' evidence>\n\n"
            "  4e. TEMPORAL TRENDS\n"
            "  ....................\n"
            "  <genome collection trend and outbreak case trend with exact slope/direction/p-value from "
            "'temporal_trend' evidence; inter-detection gap in days>\n\n"
            "  4f. INTRODUCTION / SPREAD ASSESSMENT\n"
            "  .....................................\n"
            "  Inferred origin:      <country> (probability/support: <value>)\n"
            "  Dissemination routes: <list>\n"
            "  Scenario scores:      new_introduction=<score>, re_emergence=<score>, "
            "local_persistence=<score>, with their supporting evidence strings\n"
            "  Confidence:           <level>\n\n"
            "  4g. OUTBREAK IMPACT & AFFECTED POPULATIONS\n"
            "  ...........................................\n"
            "  <For each historical outbreak caused by THIS lineage/species (sources: "
            "epidemiological_query_summary.outbreaks, comparative_outbreak_analysis, "
            "molecular_epidemiology), one line: country, year(s), cases, deaths, CFR. Then, using "
            "'demographics' and 'genomic_links', state who was affected: age group, sex, occupation/risk "
            "group, exposure setting, linked to the specific outbreak/lineage where the evidence connects "
            "them. If no demographic breakdown exists, state 'No demographic breakdown available in "
            "evidence for this lineage's outbreaks'.>\n\n"
            "  4h. VACCINES & INTERVENTIONS\n"
            "  .............................\n"
            "  Licensed/used vaccines:  <product, status, effectiveness from 'vaccines', or 'Not "
            "available in evidence'>\n"
            "  Therapeutics:            <product, status, effectiveness from 'therapeutics', or 'Not "
            "available in evidence'>\n"
            "  Other interventions:     <from 'interventions', or 'Not available in evidence'>\n\n"
            "SECTION 5: INTEGRATED EVIDENCE ASSESSMENT (EVIDENCE CONVERGENCE)\n"
            "--------------------------------------------------------------------------------\n"
            "  5a. CONVERGENT EVIDENCE — for each major conclusion (identity, evolutionary placement, "
            "phenotype significance, epidemiological context), name the independent evidence streams that "
            "agree and explain why their agreement strengthens the conclusion.\n"
            "  5b. DIVERGENT / CONFLICTING EVIDENCE — conclusions supported by only one stream, or where "
            "streams disagree; explain the conflict concretely.\n"
            "  5c. EVIDENCE STRENGTH SUMMARY — for each major conclusion, state how many independent lines "
            "of evidence support it (e.g. '3 independent streams: ...') and the resulting confidence.\n\n"
            "SECTION 6: EVIDENCE SUMMARY & KNOWLEDGE GAPS\n"
            "--------------------------------------------------------------------------------\n"
            "  6a. EVIDENCE INVENTORY (counts of literature/guidance/database sources actually used)\n"
            "  6b. CONFIDENCE SUMMARY (per analytical dimension, grounded in evidence-stated confidence)\n"
            "  6c. KNOWLEDGE GAPS & UNCERTAINTIES (explicit, categorised as scientific / surveillance / "
            "analytical)\n\n"
            "SECTION 7: GENOMIC INTELLIGENCE ASSESSMENT\n"
            "--------------------------------------------------------------------------------\n"
            "  7a. WHAT IS KNOWN  7b. WHAT HAS CHANGED vs PREVIOUS SURVEILLANCE  "
            "7c. PUBLIC HEALTH SIGNIFICANCE (biological/epidemiological significance grounded in evidence — "
            "NOT a risk tier, no MONITOR/INVESTIGATE classification, no escalation triggers).\n\n"
            "SECTION 8: VISUALISATIONS\n"
            "--------------------------------------------------------------------------------\n"
            "  <'Ref: Figure <N> – <Title>' lines per the FIGURES instructions above, each followed by the "
            "filename and a one-line grounded description; omit this section if no figures are available>\n\n"
            "================================================================================\n"
            "                         END OF REPORT\n"
            "================================================================================\n"
        )

        hard_constraints_reminder = self._hard_constraints_reminder(ascii_mode=True)

        raw = self._call(
            system_prompt=system_prompt,
            user_prompt=evidence_context_str + self._user_preface() + hard_constraints_reminder,
            temperature=0.2,
            max_tokens=10240,
            json_mode=False,
            timeout=700,
        )
        return self._enforce_ascii_contract(raw)

    @staticmethod
    def _user_preface() -> str:
        # Provide a concise, high-signal preface before the full evidence JSON so
        # smaller local models do not lose the quantitative focus in a long prompt.
        return (
            "\n\n---\n"
            "HOW TO USE THIS EVIDENCE:\n"
            "- Start by reading the 'quantitative_anchors' section. It contains the exact "
            "numbers you must use (SNP distances, percent identity, variant frequencies, "
            "lineage counts, dates, p-values, confidence scores, literature references).\n"
            "- Then read 'genomic_analyses' and 'cross_evidence_statistics' for the "
            "cross-domain links (evolutionary rate vs molecular clock, lineage history vs "
            "geography, mutation profile vs phenotype).\n"
            "- Use ONLY the evidence in this JSON. Do not supplement with training data.\n"
            "- Do not invent metric, source, or evidence-stream names; every number and fact you state "
            "must come from this JSON, written as plain prose without inline citation tags.\n"
            "\n"
            "YOUR TASK: Reproduce the EXACT plain-text template given in the system message above. "
            "Replace every '<...>' placeholder with grounded content from the evidence JSON below. "
            "Do NOT summarize the JSON, do NOT use Markdown syntax (#, **, etc.), and do NOT add sections "
            "that are not in the template. Keep the '=' and '-' divider lines exactly as shown.\n"
        )

    @staticmethod
    def _hard_constraints_reminder(ascii_mode: bool = False) -> str:
        template_rule = (
            "1. Reproduce the exact section headers and divider lines given in the template, in the exact "
            "order given, with no additional sections. Use plain text only — no Markdown (#, **, etc.).\n"
        ) if ascii_mode else (
            "1. Use EXACTLY the section headers given, in the exact order given, and no others.\n"
        )
        return (
            "\n\n---\n"
            "MANDATORY OUTPUT RULES (do not violate these):\n"
            + template_rule +
            "2. Do NOT add a 'Recommendations', 'Key Concerns', 'Suggested Actions', 'Next Steps', "
            "'Risk Classification', 'Risk Assessment', 'Response Measures', 'Escalation Triggers', "
            "'Implications', 'Evidence-Based Considerations for Decision-Making', 'Factors for the "
            "Analyst to Weigh', or any similarly action-oriented or risk-tiering section. Do not assign a "
            "risk tier (MONITOR/INVESTIGATE/HIGH PRIORITY/etc.) anywhere. Do not tell the reader what to do.\n"
            "3. Do not append inline evidence-stream citation tags like '(evidence: xxx)' or "
            "'(source: xxx)' to claims -- state facts directly as prose. Only use '(PMID: <id>)' when "
            "citing a specific literature reference, written naturally as a citation would appear in a "
            "scientific report.\n"
            "4. Do not use words like 'should', 'must', 'recommend', 'monitor for', 'strategic use of', "
            "'risk of spread', 'prepare for', or 'needs to'. Describe evidence, not next steps. "
            "Also avoid phrases such as 'further surveillance is needed', 'further research is needed', "
            "or 'more data is needed'. Instead, state exactly what evidence is missing.\n"
            "5. Do not use vague qualitative descriptors ('high', 'low', 'significant', 'notable', "
            "'concerning', 'substantial', 'major', 'important', 'widespread', 'rare', 'several', "
            "'various', 'multiple', 'primarily', 'limited', 'extensive', 'strong', 'weak') as filler; "
            "replace them with the specific number, percent, count, p-value, or metric from the "
            "evidence, stated as the sentence itself.\n"
            "6. Quote exact numbers, dates, places, and names from the evidence; do not round, "
            "generalize, or invent values.\n"
            "7. Do not add background information from your training data. Only use the evidence "
            "provided in the JSON above.\n"
            "8. Every major conclusion must explicitly name the independent evidence streams that "
            "converge to support it, or state that it rests on a single source.\n"
            "\n"
            "FINAL INSTRUCTION — REPEAT: Output MUST be the exact plain-text template from the system "
            "message, with '<...>' placeholders replaced by grounded content. No Markdown (no #, no **). "
            "No JSON summary. No extra sections. Keep the '=' and '-' divider lines exactly as shown.\n"
        )

    @staticmethod
    def _enforce_ascii_contract(text: str) -> str:
        """Guardrail for the plain-text ASCII templates (Brief/Report): strip
        any top-level section whose header matches a disallowed, risk-tiering
        or recommendation-flavored keyword, even though the prompt instructs
        the model not to add one. Local/smaller models don't always follow
        system-prompt constraints reliably, so this is a safety net on top
        of prompt engineering, not a replacement for it."""
        import re

        disallowed_keywords = (
            "risk classification", "risk assessment", "recommendation",
            "suggested action", "next step", "response measure",
            "escalation trigger", "key concern", "public health action",
            "considerations for decision-making", "factors for the analyst to weigh",
        )

        lines = text.split("\n")
        header_idx = []
        for i in range(len(lines) - 1):
            stripped = lines[i].strip()
            if (
                stripped
                and re.match(r"^-{10,}$", lines[i + 1].strip())
                and any(c.isalpha() for c in stripped)
                and stripped == stripped.upper()
            ):
                header_idx.append(i)

        if not header_idx:
            return text

        boundaries = header_idx + [len(lines)]
        keep = list(lines[: header_idx[0]])
        for idx, start in enumerate(header_idx):
            end = boundaries[idx + 1]
            header_text = lines[start].strip().lower()
            if any(kw in header_text for kw in disallowed_keywords):
                log.warning(
                    "LLM added a disallowed section ('%s'); stripping it from the output.",
                    lines[start].strip(),
                )
                continue
            keep.extend(lines[start:end])

        return "\n".join(keep).strip()

    def synthesize_genomic_intelligence(self, evidence_context_str: str) -> str:
        """Deprecated: superseded by ``synthesize_intelligence_brief`` and
        ``synthesize_intelligence_report``, which follow the project's fixed
        Brief/Report templates. Kept only for backward compatibility with
        any external callers still importing this method."""
        system_prompt = (
            "You are a senior genomic intelligence analyst. Synthesize the evidence JSON below into a "
            "narrative with sections: Executive Assessment, Genome Integrity and Analytical Confidence, "
            "Identity of the Virus, Genomic Characteristics, Evolutionary Context, Historical and "
            "Epidemiological Context, Integrated Evidence Assessment, Sources of Uncertainty, Overall "
            "Genomic Intelligence Assessment.\n\n" + self._grounding_and_reasoning_block()
        )
        raw = self._call(
            system_prompt=system_prompt,
            user_prompt=evidence_context_str + self._user_preface() + self._hard_constraints_reminder(),
            temperature=0.2,
            max_tokens=6144,
            json_mode=False,
        )
        return self._enforce_synthesis_contract(raw)

    @staticmethod
    def _enforce_synthesis_contract(text: str) -> str:
        """Guardrail: strip any recommendation/risk-flavored sections the model
        added despite instructions not to. Local/smaller models don't always
        follow system-prompt constraints reliably, so this is a safety net on
        top of prompt engineering, not a replacement for it."""
        import re

        disallowed_heading_pattern = re.compile(
            r"^#{1,6}\s*(recommendation|key concern|suggested action|next step|"
            r"risk assessment|response measure|public.health action)s?\b.*$",
            re.IGNORECASE | re.MULTILINE,
        )
        for _ in range(10):  # bounded loop: strip every disallowed section, not just the first
            match = disallowed_heading_pattern.search(text)
            if not match:
                break

            log.warning(
                "LLM added a disallowed recommendation/risk-flavored section "
                "('%s'); stripping it from the output.", match.group(0).strip()
            )
            # Cut everything from this disallowed heading to the next heading
            # of the same or shallower level, or end of text.
            heading_level = len(re.match(r"^(#{1,6})", match.group(0).lstrip()).group(1))
            next_heading_pattern = re.compile(rf"^#{{1,{heading_level}}}\s", re.MULTILINE)
            start = match.start()
            rest = text[match.end():]
            next_match = next_heading_pattern.search(rest)
            end = match.end() + next_match.start() if next_match else len(text)
            text = (text[:start] + text[end:]).strip()

        return text

    @staticmethod
    def _parse_json_response(raw: str) -> dict:
        """Parse LLM JSON response, handling common issues."""
        # Strip any markdown code fences if present
        text = raw.strip()
        if text.startswith("```"):
            # Remove ```json ... ``` or ``` ... ```
            lines = text.split("\n")
            if lines[0].startswith("```"):
                lines = lines[1:]
            if lines and lines[-1].startswith("```"):
                lines = lines[:-1]
            text = "\n".join(lines)

        try:
            return json.loads(text)
        except json.JSONDecodeError as e:
            log.error(f"LLM returned invalid JSON: {e}")
            log.debug(f"Raw response: {raw[:500]}...")
            # Return a minimal valid structure so the engine doesn't crash
            return {
                "_error": "LLM returned invalid JSON",
                "_raw_response": raw[:1000],
            }
