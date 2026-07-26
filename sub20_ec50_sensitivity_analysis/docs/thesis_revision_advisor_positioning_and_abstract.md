# Thesis Revision Notes: Advisor-Requested Positioning Changes

## Core Repositioning

The thesis should be framed as:

```text
A large-scale exploratory computational workflow for prioritizing and annotating probe-responsive cysteine-containing peptides/sites from dose-resolved chemoproteomic data.
```

It should not be framed as:

```text
A definitive target-discovery study proving direct binding mechanisms.
```

The safest thesis position is:

```text
This work prioritizes candidate probe-responsive sites and proteins, annotates their biological context, and generates structure-informed hypotheses for future validation.
```

## Advisor Revision Point 1

### Original risk

The thesis may currently sound as if it strongly concludes:

```text
These proteins are confirmed targets.
These pathways are mechanistically regulated by the probe.
Docking reveals the binding mechanism.
```

### Revised framing

Use language such as:

```text
prioritized candidate sites
annotated probe-responsive profiles
exploratory pathway-level analysis
candidate proteins for follow-up
structure-informed hypotheses
putative ligandable cysteine sites
potential transport-associated candidates
```

Avoid strong language such as:

```text
confirmed targets
validated mechanisms
directly binds
proves binding
establishes mechanism
demonstrates functional regulation
```

## Recommended Language for the Whole Thesis

### Safer verbs

Use:

```text
prioritized
annotated
suggested
highlighted
identified as candidates
was consistent with
provided support for
generated hypotheses
```

Avoid:

```text
proved
confirmed
validated
demonstrated mechanism
established direct binding
showed that the protein is a target
```

### Safer noun phrases

Use:

```text
candidate probe-responsive sites
prioritized cysteine-containing peptides
putative ligandable sites
transport-associated candidate proteins
structure-guided hypotheses
exploratory enrichment signals
```

Avoid:

```text
true targets
validated drug targets
confirmed binding sites
mechanistic targets
directly modified proteins
```

## Advisor Revision Point 2: Docking Interpretation

### Required interpretation

Docking should be described as a predictive and hypothesis-generating tool, not as experimental proof.

### Recommended thesis wording

```text
Docking and AlphaFold3-based structural inspection were used to provide structural context for selected high-confidence candidate sites. These analyses are computational predictions and should not be interpreted as direct evidence of binding, covalent modification, or functional mechanism. Rather, they were used to assess whether prioritized cysteine residues were located in structurally plausible environments for future experimental follow-up.
```

### When comparing non-covalent and covalent docking pockets

Use:

```text
The comparison between non-covalent and covalent docking poses was used to evaluate whether the predicted ligand orientation was geometrically compatible with the prioritized cysteine site.
```

Do not use:

```text
The docking result proves the ligand binds in this pocket.
```

### Safer docking conclusion

```text
The structural models suggest plausible binding-site hypotheses, but experimental validation such as competition assays, site-directed mutagenesis, or targeted MS would be required to confirm direct engagement.
```

## Advisor Revision Point 3: Abstract Must Be Especially Careful

The abstract should not overclaim target discovery or mechanism. It should emphasize:

1. large-scale computational prioritization;
2. dose-response curve modeling;
3. exploratory pathway annotation;
4. structure-guided hypothesis generation;
5. need for future validation.

## Revised Abstract Template

Use this as a safer abstract structure:

```text
Cysteine-focused chemoproteomic profiling provides a strategy for mapping ligand-responsive sites across the proteome, but large-scale dose-resolved datasets require careful computational prioritization before biological interpretation. In this thesis, I developed and applied an exploratory computational workflow to analyze TMT-based dose-resolved chemoproteomic datasets generated with the WRX-035 probe. Site-level quantitative profiles were processed using CurveCurator-based four-parameter dose-response modeling, followed by curve-quality assessment, gene-level summarization, candidate prioritization, and pathway-level annotation.

Across the 8171 dataset, 13,615 site-level profiles were successfully modeled from 15,331 raw site-level rows. Candidate prioritization was performed using curve-derived metrics including pEC50, adjusted p-value, curve fit quality, and response magnitude. Because the lowest measured concentration was 20 uM, fitted pEC50 values were interpreted within the measured concentration range, and the absence of fitted sub-20 uM EC50 estimates was treated as absence of detected below-range estimates rather than evidence that such responses were biologically absent.

Exploratory enrichment analysis suggested that the prioritized gene-level signal was distributed across transport- and RNA-localization-related biological processes rather than concentrated within a small foreground gene list. ORA did not identify significant over-represented pathways using TopN candidate lists, whereas GSEA highlighted coordinated enrichment of nucleocytoplasmic transport, intracellular protein transport, and RNA-related localization pathways. These results were interpreted as pathway-level annotations of ranked candidate genes rather than definitive evidence of pathway regulation.

Finally, AlphaFold3-based structural inspection and docking were used to provide structural context for selected candidate cysteine sites. These computational models were treated as hypothesis-generating predictions rather than direct evidence of binding or mechanism. Overall, this work establishes a reproducible exploratory workflow for prioritizing and annotating probe-responsive cysteine sites from large-scale dose-resolved chemoproteomic data and provides candidate hypotheses for future biochemical validation.
```

## Shorter Abstract Ending

If the abstract needs a shorter final sentence, use:

```text
Together, this work provides a reproducible exploratory framework for converting large-scale dose-resolved chemoproteomic measurements into prioritized candidate sites and structure-informed hypotheses for future experimental validation.
```

## Revised Thesis-Wide Claim Hierarchy

### Strongest claim allowed

```text
The workflow prioritized and annotated candidate probe-responsive cysteine sites.
```

### Moderate claim allowed

```text
Pathway analysis suggested a distributed enrichment pattern involving RNA localization, nucleocytoplasmic transport, and intracellular protein transport.
```

### Structural claim allowed

```text
Docking and structural modeling provided plausible hypotheses for selected candidate sites.
```

### Claims to avoid

```text
WRX-035 directly targets IPO5/IPO7/PSME1/DHX57.
The docking models confirm covalent binding.
The enriched pathways are mechanistically regulated by WRX-035.
The identified proteins are validated targets.
```

## Revised Results Tone

Instead of:

```text
The analysis identified WRX-035 targets involved in nuclear transport.
```

Use:

```text
The analysis prioritized candidate probe-responsive proteins whose gene-level rankings were enriched for transport-related annotations.
```

Instead of:

```text
Docking confirmed that WRX-035 binds to IPO5 Cys682.
```

Use:

```text
Docking generated a structurally plausible pose near IPO5 Cys682, supporting this site as a candidate for future experimental validation.
```

Instead of:

```text
GSEA revealed the mechanism of WRX-035.
```

Use:

```text
GSEA provided pathway-level annotation suggesting that transport- and RNA-localization-related genes were broadly shifted in the ranked candidate list.
```

## Revised Discussion Framing

Use:

```text
The present analysis should be interpreted as an exploratory prioritization workflow. The combination of curve fitting, enrichment analysis, and structural modeling narrowed a large chemoproteomic dataset into a smaller set of annotated candidate sites and proteins. However, the results do not establish direct binding or functional mechanism. Experimental follow-up, including orthogonal competition assays, targeted mass spectrometry, mutagenesis of candidate cysteine residues, and biochemical validation, would be required to confirm direct probe engagement and functional relevance.
```

## Integration with 20 uM Boundary Limitation

The 20 uM boundary issue reinforces the same conservative thesis position:

```text
pEC50 values were used primarily as exploratory ranking metrics rather than definitive potency measurements for all profiles.
```

Recommended combined limitation sentence:

```text
Because the experiment lacked an untreated channel and concentrations below 20 uM, and because docking was used only as a computational prediction, the resulting candidate sites should be interpreted as prioritized hypotheses rather than validated targets or mechanisms.
```

## Final Thesis Position

The final thesis should consistently communicate:

```text
This project developed and applied a reproducible exploratory workflow for large-scale prioritization and annotation of candidate probe-responsive cysteine sites. The workflow integrates dose-response modeling, enrichment analysis, and structural prediction to generate biologically interpretable hypotheses, but experimental validation is required before drawing firm conclusions about direct binding, target identity, or mechanism.
```
