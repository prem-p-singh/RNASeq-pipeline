# Integration validation — 2026-10-03

This records orchestration tests, not a new blanket scientific qualification or a published 3.0 release.

## Source and fixes

- Base: release-3.0 commit `5936319a355f12f1418c3117ae8bc5a8c29b2e10`.
- Deployed scientific-source fingerprint, including the recorded bootstrap patch: `c0927cf7f31e0fc603c39eec9b432f9efcdcbbc53971d422dfe1ba2b2de6750c`.
- The bridge now preserves an already-active Conda interpreter ahead of base Conda. The first Linux attempt exposed this issue and was not a pass.
- `patches/bootstrap-path.patch` explicitly puts the environment being verified first on PATH. The second attempt exposed STAR's samtools shadowing the UMI module during verification; that failed run is not counted as a pass.

## Passed checks

| Check | Evidence and boundary |
|---|---|
| Import | Corrected workflow imported into n8n 2.25.7, both isolated test database and the user's existing project. Original workflow retained. |
| Actual n8n control flow | CLI execution in the isolated database reached `Review results`, status `success`, exit 0. SSH nodes were replaced with simulated responses and the wait reduced to one second. This tests branching/polling, **not SSH authentication**. |
| Bridge lifecycle | Local and Linux checks: completion/failure, same-ID reconnect, changed-request refusal, project locking, interruption, source fingerprint and project identity checks. The fixture now waits for worker exit before the next case to avoid racing final lock release. |
| Code nodes | Settings, source gate, response identity/error handling and route presentation passed with Node.js. |
| Intake | Release `tests/check_intake.py` passed in Linux job 39402586: current workbook schema, bulk/TAGseq/sRNA setup, previous schema handling, advanced choices, revision reporting and explicit apply. These are the existing fixture tests, not a native Excel UI session. |
| Independent raw QC | Job 39402586 ran real fastp through the bridge; successful report and unchanged original reads. This run predates the bootstrap patch; the core runtime path was exercised. |
| QuantSeq FWD and REV | Real launcher and Snakemake dry runs passed in job 39402616. |
| QuantSeq FWD-UMI | Dry run and real synthetic analysis passed through the bridge in job 39402616; gene A has 20 molecules in each of six samples, gene B zero. |
| Animal small RNA | Dry run and real synthetic analysis passed through the bridge in job 39402616; miR-1 counts 5/3/6 and miR-3 counts 2/4/0 match the fixture. |
| Plant small RNA | Real launcher/DAG dry run passed in job 39402616. No new plant scientific analysis was run by this integration check. |
| Bulk Salmon | Real launcher/DAG dry run passed in job 39402616. |
| Bulk STAR | Focused real launcher/DAG dry run passed in job 39403288, completed exit 0. |

Job 39402616 did **not** pass as a whole: its final bulk-STAR fixture placed expected library type in a sample-sheet column rather than `samples.expected_libtype`. Preflight correctly rejected it. The fixture was corrected without weakening pipeline validation. Focused job 39403288 passed after that correction. Together, the recorded checks cover all seven route/kit combinations; this is not represented as a single clean full-suite run.

## Live OpenSSH connection verified

The user authorized using the existing `ssh farm` alias. The running native n8n instance already enabled Execute Command (no instance security setting was changed). The OpenSSH variant was executed by n8n 2.25.7 in an isolated CLI database against the actual analysis host, without simulated SSH responses. It verified the deployed source, retrieved the completed bulk-STAR run, reached `Review results`, reported `star_counts` / `succeeded`, and finished with n8n status `success` and exit 0.

The identical host-configured workflow was imported into the user's existing n8n project as **RNASeq 2.1 + 3.0 — FARM via OpenSSH**, retaining status-only mode. No private key was copied into n8n, no SSH credential was created, and host-key checks remain enabled. The portable OpenSSH export uses a configurable `analysis-host` alias instead of assuming a cluster name. The credential-based export remains available for other deployments.

This live test reconnects to completed work; it does not claim a new scientific dataset was launched from the browser. Launch, polling and scientific execution were exercised separately in the Linux bridge tests above. For a new study, configure its host paths and appropriate executor before changing the action to start.

Scientific route tests above use `executor=local` within a SLURM worker allocation. A new distributed `executor=slurm` n8n submission was not exercised here. Existing scheduler behavior is delegated to `submit.sh` and its configured profiles. The full release scientific suite was not rerun by this integration task.

Detailed local evidence is retained in ignored `docs/evidence/n8n-20261003/`; it is not included in the public source archive.

## Parent integration gate — 2026-10-05

The parent remains schema 4; this handoff does not implement the schema-7 routes
in that parent. `tests/run_all.sh --n8n` explicitly dispatches the node, bridge,
and reviewed-source release checks. Missing release prerequisites are SKIP, not
PASS, and prevent `--strict` qualification.

The audit's local interpreter failure was reproduced. The agent runtime supplied
`BASH_ENV` with a fixed PATH export, overriding the bridge's active Conda prefix
when its fixture launcher started Bash. The unchanged check passed with that
hook absent. The isolated fixture now clears `BASH_ENV` in its own test process;
the production bridge and the exact interpreter-path assertion are unchanged.
This is a test-environment correction, not evidence of a new Linux release run.

Measured local gate result: **2 passed, 0 failed, 1 skipped**, developer exit 0;
strict mode returns 1 with `RELEASE BLOCKED`. The release smoke skipped because
this host is macOS and has no configured reviewed deployment or locked Linux
module prefixes. Bash/Python syntax checks passed. A simulated prerequisite
check rejected the schema-4 parent before launching anything. The rebuilt bundle
reproduced the recorded scientific fingerprint and contained no Python caches.
No new Linux smoke, full scientific suite, remote submission, or release tag was
performed for this change.
