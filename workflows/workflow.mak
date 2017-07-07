
# Input parameters
# ================
#
path_nimbus      := $(shell echo $$PATH_NIMBUS)
path_samtools    := $(shell echo $$PATH_SAMTOOLS)
path_python      := $(shell echo $$PATH_PYTHON)
path_zcat        := $(shell echo $$PATH_ZCAT)
path_gzip        := $(shell echo $$PATH_GZIP)
path_annovar     := $(shell echo $$PATH_ANNOVAR)
annovar_db       := $(shell echo $$ANNOVARDB)
amplicon_design  := $(shell echo $$AMPLICON_DESIGN)
genome_reference := $(shell echo $$GENOME_REFERENCE)

# Minor parameters
# -------------------

# Trim
minimum_matches      = 2
maximum_mismatches   = 2
minimum_readlength   = 50
adapter_trim_options = --maximum-mismatches $(maximum_mismatches) --minimum-matches $(minimum_matches) --minimum-bases-remaining $(minimum_readlength)
adapters             = -a AGATCGGAAGAG -a CTGTCTCTTATA

# Alignment
keysize = 6
workers = 4

# SAMtools
sort_options   := --threads 4
library_label  := Amplicon
platform_label := ILLUMINA
center_label   := Unknown

# Filter
maxscore := 3

# Calling
min_mapping_qual  := 0
max_reads         := 1000000
max_alleles       := 65536
min_allelic_depth := 1
non_reference     := 1

# Var
minimum_quality_frequency := 0.2
minimum_quality := 10

# Annovar options
version   := hg19
protocol  := refGene,genomicSuperDups,avsnp147,gerp++gt2
operation := g,r,f,f

# Define the target files
# =======================
#
fileinput     := $(wildcard *_R1.fastq) $(wildcard *_R1.fastq.gz) $(wildcard *_R1_001.fastq) $(wildcard *_R1_001.fastq.gz) $(wildcard *.srt.bam) 
filebase      := $(sort $(foreach entry, $(fileinput), $(shell echo $(entry) | sed 's/_R1_001.fastq.*$$//'| sed 's/_R1.fastq.*$$//' | sed 's/.srt.bam.*$$//')))
bamfiles      := $(patsubst %, %.srt.bam, $(filebase))
bamflags      := $(patsubst %, %.flagstat.txt, $(filebase))
bampassed     := $(patsubst %, %.passed.bam, $(filebase))
bamdiscarded  := $(patsubst %, %.discarded.bam, $(filebase))
blckfiles     := $(patsubst %, %.srt.blck, $(filebase))
blckpassed    := $(patsubst %, %.passed.blck, $(filebase))
blckdiscarded := $(patsubst %, %.discarded.blck, $(filebase))

# Control flow
# ============
#
all: alignment flagstats filter coverage calling annotation visualization table

alignment: $(bamfiles)

coverage: $(blckfiles) $(blckpassed) $(blckdiscarded)

flagstats: $(bamflags)

filter: $(bampassed) $(bamdiscarded)

calling: combined.var.gz combined.txt

annotation: combined.hg19_multianno.txt

visualization: combined.mut.txt

table: combined.hg19.anno.txt

# Recipes
# =======

# Trim files
# ----------
%.tr.fastq: %_001.fastq.gz
	mkdir -p logs
	$(path_zcat)/zcat $*_001.fastq.gz | $(path_nimbus)/bin/nimbus_trim \
		-i - \
		$(adapter_trim_options) \
		$(adapters) \
		-o $*.tr.fastq 2>> logs/$*.trim.errors.log >> logs/$*.messages.errors.log

%.tr.fastq: %.fastq.gz
	mkdir -p logs
	$(path_zcat)/zcat $*.fastq.gz | $(path_nimbus)/bin/nimbus_trim \
		-i - \
		$(adapter_trim_options) \
		$(adapters) \
		-o $*.tr.fastq 2>> logs/$*.trim.errors.log >> logs/$*.messages.errors.log

%.tr.fastq: %_001.fastq
	mkdir -p logs
	$(path_nimbus)/bin/nimbus_trim \
		-i $*_001.fastq \
		$(adapter_trim_options) \
		$(adapters) \
		-o $*.tr.fastq 2>> logs/$*.trim.errors.log >> logs/$*.messages.errors.log

%.tr.fastq: %.fastq
	mkdir -p logs
	$(path_nimbus)/bin/nimbus_trim \
		-i $*.fastq \
		$(adapter_trim_options) \
		$(adapters) \
		-o $*.tr.fastq 2>> logs/$*.trim.errors.log >> logs/$*.messages.errors.log

# Nimbus alignment
# ----------------
%.sam: $(amplicon_design) $(genome_reference) %_R1.tr.fastq %_R2.tr.fastq
	mkdir -p logs
	$(path_nimbus)/bin/nimbus_align align \
		--forward $*_R1.tr.fastq \
		--reverse $*_R2.tr.fastq \
		--design $(amplicon_design) \
		--fasta $(genome_reference) \
		--key-size $(keysize) \
		--workers $(workers) \
		--maximum-amplicons 1000 \
		--sam $*.sam 2>> logs/$*.nimbus.errors.log >> logs/$*.nimbus.messages.log

# SAMtools processing
# -------------------
%.srt.bam: %.sam
	sample=$$(echo $* | sed 's/_.*$$//') ; \
	path=$$(readlink -f $*.srt.bam) ; \
	$(path_samtools)/samtools addreplacerg \
		-r "ID:$${sample}" \
		-r "CN:$(center_label)" \
		-r "LB:$(library_label)" \
		-r "SM:$${sample}" \
		-r "PL:$(platform_label)" \
		-r "DS:$${path}" \
		$*.sam | $(path_samtools)/samtools sort \
			$(sort_options) \
			-T $*_tmp_sort \
			-o $*.srt.bam \
			- ;

%.flagstat.txt: %.srt.bam
	$(path_samtools)/samtools flagstat $*.srt.bam > $*.flagstat.txt

# Filter bad reads
# ----------------
%.temp.bam %.discarded.bam: %.srt.bam
	mkdir -p logs
	$(path_python)/python $(path_nimbus)/scripts/nimbus_filter.py \
		--score $(maxscore) \
		--input $*.srt.bam \
		--output $*.temp.bam \
		--discarded $*.discarded.bam 2>> logs/$*.filter.errors.log >> logs/$*.filter.messages.log

%.passed.bam: %.temp.bam
	$(path_samtools)/samtools sort \
		$(sort_options) \
		-T $*_tmp_passed_sort \
		-o $*.passed.bam \
		$*.temp.bam

# Count the reads per amplicon
# ----------------------------
%.blck: $(amplicon_design) %.bam
	mkdir -p logs
	$(path_python)/python $(path_nimbus)/scripts/nimbus_count.py \
		--input $*.bam \
		--design $(amplicon_design) \
		--output $*.blck >> logs/$*.count.messages.log 2>> logs/$*.count.errors.log

# Call variants
# -------------
combined.var.gz: $(genome_reference) $(bampassed)
	mkdir -p logs
	$(path_nimbus)/bin/nimbus_call \
		--bam=`echo $(bampassed) | sed 's/ / --bam=/g'` \
		--fasta=$(genome_reference) \
		--info=am \
		--minimum-mapping-quality=$(min_mapping_qual) \
		--maximum-number-of-reads-in-pileup=$(max_reads) \
		--maximum-number-of-alleles=$(max_alleles) \
		--minimum-allelic-depth=$(min_allelic_depth) \
		--report-only-non-reference-alleles=$(non_reference) | \
			$(path_gzip)/gzip -c > combined.var.gz 2>> logs/combined.call.errors.log

combined.txt: combined.var.gz
	mkdir -p logs
	 ${path_zcat}/zcat combined.var.gz | \
		$(path_python)/python $(path_nimbus)/scripts/nimbus_var.py \
			--input - \
			--remove-N \
			--minimum-quality-frequency $(minimum_quality_frequency) \
			--minimum-quality $(minimum_quality) \
			--output combined.txt 2>> logs/combined.var.errors.log

# Annovar
# -------
combined.hg19_multianno.txt: $(annovar_db) combined.txt
	mkdir -p logs
	mkdir -p combined_annovar_tmp
	$(path_annovar)/table_annovar.pl \
		-buildver $(version) \
		-protocol $(protocol) \
		-operation $(operation) \
		-nastring - \
		--otherinfo \
		--outfile combined \
		--tempdir combined_annovar_tmp \
		combined.txt \
		$(annovar_db) 2>> logs/combined.annovar.errors.log >> logs/combined.annovar.messages.log

# Make mutation file
# ------------------
combined.mut.txt: combined.hg19_multianno.txt
	mkdir -p logs
	$(path_python)/python $(path_nimbus)/scripts/nimbus_mut.py \
		--filter 'Quality_frequency:<$(minimum_quality_frequency),->' \
		--input combined.hg19_multianno.txt \
		--output combined.mut.txt 2>> logs/combined.mut.errors.log >> logs/combined.mut.messages.log

combined.hg19.anno.txt: combined.hg19_multianno.txt
	mkdir -p logs
	$(path_python)/python $(path_nimbus)/scripts/nimbus_table.py \
		--input combined.hg19_multianno.txt \
		--output combined.hg19.anno.txt 2>> logs/combined.table.errors.log >> logs/combined.table.messages.log
