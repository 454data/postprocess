# postprocess

## Installation
Set the path-constants in the script header to point at:
```
USEARCH    usearch v5.2.32
HMMSEARCH  hmmsearch v3.1
BLASTN     blastn v2.2.25
RESOURCES  contents of resources.zip
```

## Usage
```
Usage: postprocess_basecalls.py [options] <basecalls_file>

<basecalls_file> File with '.basecalls'-extension produced by Multipass.

This script post-processes amplicon sequenced DBLa-tags basecalled using
Multipass.

Options:
  --version      show program's version number and exit
  -h, --help     show this help message and exit
  -v             print verbose information [False]
  -R RESULT_DIR  directory for result files [<basecall_file>.postprocess]
  -i             remove 3D7 sequences using BLASTN [False]
  -m METHOD      basecalling method. 0=Multipass, 1=Multipass_FRF  [1]
```

## Result files
```
## Basecalls ##
.bc.fas
## Clustering ##
.clu.fas
.clu.uc
.clus.fas
## Chimeras de novo ##
.clus.uchimealns
.clus.uchimeout
.clus.ch.fas
.clus.nc.fas
## Chimeras database mode ##
.clus.nc.db
.clus.uchimealnsself
.clus.uchimeoutself
.clus.chself.fas
.clus.ncself.fas
## Coverage trim ##
.clust.fas
## Non-DBLa ##
.atag.hmmsearch
.nonatag.fas
## DBLb-tags ##
.btag.hmmsearch
.btag.fas
## 3D7 removed ##
.3d7atags.fas
.3d7nonatags.fas
## Final cleaned DBLa-tags ##
.clean.fas
```
