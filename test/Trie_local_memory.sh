#!/bin/bash
# Trie script

keyFile=300k_rule_subset-evict-noise-withID-

ipfileFolder=~/Research/PCAP/mapped_caida_trace/

memSize=(100 200 250 300 320 340 360 380 410 450)

feedbackPortion=(0.01 0.05 0.1 0.2 0.4 0.6)

interval=(2000000 5000000 10000000 20000000 50000000 100000000)

blackKeySize=(500 1000 1500 2000 4000 8000)

for i in $(seq 0 0)

do
	export OMP_NUM_THREADS=8

	./trieNoiseMain ${keyFile} ${memSize[i]} ${ipfileFolder} ${feedbackPortion[2]} ${interval[0]} ${blackKeySize[1]}

done






