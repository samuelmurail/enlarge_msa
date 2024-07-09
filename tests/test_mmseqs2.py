#!/usr/bin/env python

"""Tests for `enlarge_msa` package."""

import pytest
import sys
sys.path.insert(0, "src/")


from enlarge_msa import mmseqs2


def test_mmseqs2_create_full_alignement(tmp_path):
    
        
    original_seq = "GPEKEWVEQDEPGVYITLTALAGGARDLKRVRFSRKRFSEIQAEQWWADNRGRVYEQYNVRMV"
    original_seq += ":GPKWVKTDSDFIVLEI"
    original_seq_list = original_seq.split(":")

    new_seq = "EIVKLRTQVDSLTQKCEFQEVELQNSVKKTQEALALAEEESAKSRAAKEAIKSLIAQLKDVAEKLPPGES"\
              "VKLACLQNGLDQNGFHFPEENGFHPSRSESMTSSISSVAPFDFAFANASWSNLQSPKQTPRASERNSNAY"\
              "PADPRLSSSGSVISERIEPFQFQNNSDNGSSQTGVNNTNGPEKEWVEQDEPGVYITLTALAGGARDLKRV"\
              "RFSRKRFSEIQAEQWWADNRGRVYEQYNVRMVEKSTASQTHRDRDEEEEDIPH"
    new_seq += ":DEDIEKTLSVILDKCKDICAEKSKKSIGKKSISFLLKKMFVCRSGFAPTPSLRDTLQESRMEKLLRTMLH"\
              "KKLYTQNNSRAPVLKKCLENKKSIKKRNEDEAEERIDEGPKWVKTDSDFIVLEI"

    new_seq_list = new_seq.split(":")

    mmseqs2.create_full_alignement(
            '6LV0',
            original_seq_list,
            new_seq_list,
            out_dir=tmp_path,
            #out_dir='tmp',
            )
    
    with open(tmp_path.joinpath("6LV0.a3m"), "r") as f:
        result = f.read()
    
    seq_num = 0
    seq_num_type = {}
    seq_type = ""
    for line in result.split("\n"):
        print(f'#{line}#')
        if line.startswith(">"):
            seq_num += 1
            if line[1:] == "101\t102":
                seq_type = "paired"
                # print("paired")
            if len(line[1:]) == 3:
                # print("unpaired", int(line[1:-1]))
                seq_type = line[1:]
            if seq_type not in seq_num_type:
                seq_num_type[seq_type] = 0
            seq_num_type[seq_type] += 1

    print(seq_num_type)
    assert seq_num_type['paired'] == 30
    assert    seq_num_type['101'] == 2171
    assert    seq_num_type['102'] == 250

