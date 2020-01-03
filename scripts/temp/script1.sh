
#!/bin/bash
ls fa | xargs -P40 -I {} Rscript CRISP1.R {}
