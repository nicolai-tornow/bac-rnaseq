# Reference bundles

Committed inputs (Rock-lab annotation):

- `mabs/NC_010397.1.{fasta,gff3}` — *Mycobacterium abscessus* ATCC19977 chromosome.
- `mtb/H37RvBD.{fasta,gff3}` — *M. tuberculosis* H37Rv. NOTE: the GFF seqid is
  `H37RvBD` but the FASTA header is `NC_018143.2`; `build-refs` remaps the GFF
  seqid to the FASTA one so the SAF `Chr` matches the alignment references.

Built at setup (git-ignored): `<species>/index/ref.*`, `<fasta>.fai`, `labels.saf`.

Custom genomes: the user supplies FASTA + GFF3 at run time; nothing is committed here.
