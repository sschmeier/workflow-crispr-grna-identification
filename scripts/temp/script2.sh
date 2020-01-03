#!/bin/bash

ls script1out | xargs -P40 -I {} Rscript CRISP2.R {}
