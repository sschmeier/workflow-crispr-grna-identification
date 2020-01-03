#!/bin/bash

ls script2out/top25fa | xargs -P83 -I {} Rscript CRISP3.R {}
