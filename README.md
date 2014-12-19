reference_corrector
===================

Just a small step in reference correction for little bioinformaticians, but a giant leap for a mankind. This source code uses deBruijn graph and other data structures from SPAdes (>=3.1.1) source code.

Usage:

  1. download latest version of SPAdes

  2. add or replace files:
    - from /debruijn to src/debruijn/ folder
    - from /configs to configs/ folder
  
  3. build SPAdes as usual.

Code from /debruijn folder adds new stage to SPAdes and tries to correct reference.
