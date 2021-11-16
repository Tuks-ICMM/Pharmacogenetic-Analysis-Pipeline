# Pipeline Configuration:

Welcome to the configuration folder. This folder conatins all the configuration ifnormation that controlls our understanding of the environment you are running thsi pipeline on, as well as default scheduler profiles for several common Batch Scheduler environments.

> We currently only support PBS-Torque out the box. You are more than welcome to add in support for your own. If you do, we ask if you would be willing to PR the new profiles to contribute to this repository for the benefit of others.

## Format spesification
A `config.schema.json` file has been created to validate yoru configuration file to ensure it is properly formatted. **You do not need to go through this yourself**. The pipeline will automatically perform validation for you and inform you if the `config.json` file is not formatted correctly.

Here is a sample, using the genes CYP2A6, CYP2B6 and UGT2B7 with the 1000 Genomes dataset:

```json
{
  "reference-genomes": [
    {
      "version": "GRCh38",
      "file": "GRCh38.fa.gz"
    }
  ],
  "samples": [
    {
      "name": "1000g",
      "reference_version": "GRCh37",
      "exclusions": [
        "rs199808813"
      ]
    }
  ],
  "locations": [
    {
      "name": "CYP2A6",
      "transcripts": [
        "NM_000762.6",
        "ENST00000600495.1",
        "ENST00000596719.5",
        "ENST00000599960.1"
      ],
      "chromosome": 19,
      "strand": "-1",
      "start": 40842850,
      "stop": 40851138
    },
    {
      "name": "CYP2B6",
      "transcripts": [
        "NM_000767.5",
        "ENST00000593831.1",
        "ENST00000598834.2",
        "ENST00000597612.1",
        "ENST00000594187.1"
      ],
      "chromosome": 19,
      "strand": "1",
      "start": 40988570,
      "stop": 41021110
    },
    {
      "name": "UGT2B7",
      "transcripts": [
        "NM_001074.4",
        "ENST00000508661.5",
        "ENST00000622664.1",
        "ENST00000502942.5",
        "ENST00000509763.1"
      ],
      "chromosome": 4,
      "strand": "1",
      "start": 69045214,
      "stop": 69112987
    }
  ],
  "cluster": {
    "clusters": [
      "SUPER",
      "SUB"
    ],
    "file": "Clusters.xlsx"
  },
  "environment": {
    "email": {
      "address": "graeme.ford@tuks.co.za",
      "conditions": [
        "o",
        "e"
      ]
    },
    "working-directory": "/nlustre/users/graeme/pipeline-2020",
    "queues": [
      {
        "name": "long",
        "walltime": "900:00:00",
        "memory": "128G"
      }
    ]
  }
}
```

## Notes
There are a few things to take note of as indicated below.

> More comprehensive descriptions comming soon