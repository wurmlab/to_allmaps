### Map optical maps

``
# Run hybrid scaffolding in IrysView and copy over the output from prometheus.
scp -r prometheus:/scratch/local/estolle/irys/data/2019-06-06-priyam/Sinvicta_SB/hybridscaffold/output/hybrid_scaffolds results/2019-06-11-hybrid_scaffolds

# Link the XMAP file obtained from hybrid scaffolding step to tmp/
ln -s $(readlink -f results/*/*_NGScontigs_*.xmap) tmp/BNG.xmap

# Convert XMAP to BED and filter.
./scripts/make_optical_map.sh
```
