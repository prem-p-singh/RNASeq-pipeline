# Scheduling and storage

The scientific route follows the declared assay, requested objectives, available implementations and dependence structure. Sample count does not switch the DE model or invent an unsupported measurement route. See [recommendation rules](config/recommendation_rules.yaml).

The launcher currently chooses a small/medium/large FARM profile for scheduler defaults. These are scheduling conveniences, not scientific categories. Each rule requests its own CPU, memory and time allocation. `scripts/plan_resources.py` derives a storage-compatible upper bound on concurrency; `hpc.samples_in_flight` can lower it further.

There is no fixed platform storage cap. Planning groups paths on the same filesystem and accounts for new retained output, retained remote downloads, missing references, index workspace, environments and concurrent temporary files. Existing local reads and detected caches are not charged as new allocations. Known input bytes influence the temporary estimate; remote sizes that cannot be measured use labelled assumptions.

The planner compares a conservative peak estimate plus margin/reserve with free space and any explicit remaining quota. It lowers concurrency first. If one job or the final retained output cannot fit, it reports the shortfall and stops. It does not silently change the scientific analysis or delete original FASTQs to fit.

The shipped coefficients are still assumptions. There is no continuous mid-run capacity monitor, remote HEAD-based sizing, or calibrated large-cohort memory model. Account quota discovery is site-dependent; `df` alone can overstate user capacity on shared storage. Supply a real remaining quota when required and select project/cache locations with enough capacity.

Use [README.md](README.md) for environment setup and launch commands, and [docs/RELEASE.md](docs/RELEASE.md) for validation evidence and limitations.
