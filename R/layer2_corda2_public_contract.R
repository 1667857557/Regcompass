# Compatibility shim.
#
# CORDA2 parameter parsing, cache construction, output metadata and persistence
# are implemented directly in `layer2_corda_runtime.R`. Worker-pool ownership is
# implemented in `layer2_corda_pool_lifecycle.R`. This file intentionally adds
# no function override; it remains in Collate only so development installations
# created from earlier PR #254 revisions keep a stable source-file list.
