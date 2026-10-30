#include "common/acsConfig.hpp"

map<string, string> ACSConfig::docs = {

    {"product_basis_mode", R"config(Explicit product-basis model for zhang_full_rank.

SINGLE_S_BASIS_NATIVE is the default ONLY in the isolated P0/P1 candidate source.
The current committed estimator S-basis is the sole basis authority; no independent
product tree, core reselection, or legacy coordinate fallback is permitted.
Selecting this mode is NOT an integer-fix decision or an AR publication certificate.

The old product-tree path is available only through the exact explicit value
LEGACY_PRODUCT_TREE. Unknown, empty, lower-case, whitespace-padded and aliased values
are rejected. An omitted value emits a candidate-default migration warning.

Old run files and checkpoints are not migrated by this default change. A checkpoint
must carry and match the product-basis mode/schema; a legacy or unmarked checkpoint
must not be restored as a native single-S-basis state.
)config"},

    {"outputs", R"config(Specifies options to enable outputs and specify file locations.

Each section typically contains an option to `output` the filetype, and a `directory` to place the files named `filename`, along with any ancillary options.
)config"},

    {"inputs",
     R"config(This section of the yaml file specifies the lists of files to be used for general metadata inputs, and inputs of external product data.
)config"},

    {"gnss_observations",
     R"config(This section specifies the sources of observation data to be used in positioning.


There are numerous ways that the `pea` can access GNSS observations to process.
You can specify individual files to process, set it up so that it will search a particular directory, or you can use a command line flag `--rnx <rnxfilename>` to add an additional file to process.

The data should be uncompressed rinex (gunzipped, and not in hatanaka format), or RTCM3 formatted binary data.


It may consist of RINEX files, or RTCM streams or files, which are specified as follows:

```
	gnss_observations:
		root_stations_directory: /data/acs/ginan/examples/data
		rnx_inputs:
			- "ALIC*.rnx"
			- "BAKO*.rnx"

		#rtcm_inputs:
		#	- "*-OBS.rtcm3"

		#streams:
		# - "https://<USER>:<PASS>@ntrip.data.gnss.ga.gov.au/ALIC00AUS0"
```

The first 4 characters of the filename are used as the receiver ID.

If multiple files are supplied with the same ID, they are all processed in sequence - according to the epoch times specified within the files. In this case, it is advisible to correctly specify the start_epoch for the filter, or the first epoch in the first file will likely be used.

)config"},

    {"satellite_data",
     R"config(This section specifies sources of ephemerides and other satellite data.)config"},

    {"trace", "Trace files are used to document processing"},

    {"ssr", "Values derived from applying received corrections to broadcast ephemeris"},
    {"broadcast", "Values derived from broadcast ephemeris streams/files"},
    {"precise", "Values derived from file-based products such as SP3/CLK/OBX"},
    {"kalman", "Values estimated internally by the kalman filter"},

    {"cost",
     "COST format files are used to export troposhere products, such as ZTD and delay gradients."},
    {"trop_sinex",
     "Troposphere SINEX files are used to export troposhere products, such as ZTD and delay "
     "gradients."},
    {"slr_obs",
     "SLR_OBS files are used as temporary files to arrange SLR observations by time. SLR "
     "observations are taken from "
     "CRD files, which are not strictly in time-order)."},
    {"slr_options",
     "This section controls how Satellite Laser Ranging (SLR) observations are handled."},
    {"ssr_outputs",
     "This section specifies how State State Representation (SSR) corrections are calculated "
     "before being published to "
     "an NTRIP caster."},
    {"rtcm_inputs",
     "This section specifies how State State Representation (SSR) corrections are applied after "
     "they are downloaded "
     "from an NTRIP caster."},

};
