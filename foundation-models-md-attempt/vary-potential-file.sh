while IFS= read -r path; do
    base="$(basename "$path" .pt)"
    sed -e "s|POTENTIAL_FILEE|$path|g" -e "s|POTENTIAL_FILE_BASE|$base|g" input-lammps-144atoms-min-only.in > "input-$base.prod.in"
    mpirun -np 1 lmp -sf kk -k on g 1 -in input-$base.prod.in | tee log-$base.log
done < /home/cganley2/scitech2026/mace-models/foundation-models/lammps-potentials/foundation-models-paths.txt