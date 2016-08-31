{
  "class": "CommandLineTool",
  "outputs": [
    {
      "type": [
        "null",
        "File"
      ],
      "outputBinding": {
        "sbg:metadata": {
          "file_type": "text"
        },
        "secondaryFiles": [
          "^.tsv"
        ],
        "loadContents": false,
        "glob": {
          "class": "Expression",
          "script": "\"./targqc_results/summary.b64html\"",
          "engine": "#cwl-js-engine"
        }
      },
      "label": "Summary report",
      "description": "A table with sample-wide statistics for each sample",
      "sbg:fileTypes": "BASE64HTML",
      "id": "#summary_report"
    },
    {
      "type": [
        "null",
        "File"
      ],
      "outputBinding": {
        "loadContents": true,
        "glob": {
          "class": "Expression",
          "script": "\"./targqc_results/regions.tsv\"",
          "engine": "#cwl-js-engine"
        }
      },
      "label": "Region-based report",
      "description": "Statistics for each input region, as well as overlapping exon (or for the full exome, if no BED is provided).",
      "sbg:fileTypes": "TSV",
      "id": "#region_report"
    }
  ],
  "stdout": "",
  "description": "",
  "sbg:latestRevision": 24,
  "sbg:contributors": [
    "vladsaveliev"
  ],
  "sbg:project": "vladsaveliev/testproject",
  "temporaryFailCodes": [],
  "sbg:createdBy": "vladsaveliev",
  "hints": [
    {
      "class": "sbg:CPURequirement",
      "value": {
        "class": "Expression",
        "script": "$job.inputs.threads ? $job.inputs.threads : 1",
        "engine": "#cwl-js-engine"
      }
    },
    {
      "class": "sbg:MemRequirement",
      "value": 1000
    },
    {
      "dockerImageId": "",
      "class": "DockerRequirement",
      "dockerPull": "images.sbgenomics.com/vladsaveliev/targqc:v1"
    }
  ],
  "sbg:createdOn": 1472150265,
  "baseCommand": [
    "targqc"
  ],
  "sbg:image_url": null,
  "sbg:modifiedBy": "vladsaveliev",
  "sbg:revisionsInfo": [
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 0,
      "sbg:modifiedOn": 1472150265,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 1,
      "sbg:modifiedOn": 1472150735,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 2,
      "sbg:modifiedOn": 1472150809,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 3,
      "sbg:modifiedOn": 1472151011,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 4,
      "sbg:modifiedOn": 1472151124,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 5,
      "sbg:modifiedOn": 1472498552,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 6,
      "sbg:modifiedOn": 1472499819,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 7,
      "sbg:modifiedOn": 1472499857,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": "Add sbg docker hub url",
      "sbg:revision": 8,
      "sbg:modifiedOn": 1472504973,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 9,
      "sbg:modifiedOn": 1472507513,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 10,
      "sbg:modifiedOn": 1472507522,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 11,
      "sbg:modifiedOn": 1472510326,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 12,
      "sbg:modifiedOn": 1472516669,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 13,
      "sbg:modifiedOn": 1472516760,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 14,
      "sbg:modifiedOn": 1472516952,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 15,
      "sbg:modifiedOn": 1472551222,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 16,
      "sbg:modifiedOn": 1472551256,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 17,
      "sbg:modifiedOn": 1472552981,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 18,
      "sbg:modifiedOn": 1472553693,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 19,
      "sbg:modifiedOn": 1472554694,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 20,
      "sbg:modifiedOn": 1472555449,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 21,
      "sbg:modifiedOn": 1472556172,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": "Add html_to_b64html.py",
      "sbg:revision": 22,
      "sbg:modifiedOn": 1472565868,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 23,
      "sbg:modifiedOn": 1472566609,
      "sbg:modifiedBy": "vladsaveliev"
    },
    {
      "sbg:revisionNotes": null,
      "sbg:revision": 24,
      "sbg:modifiedOn": 1472567627,
      "sbg:modifiedBy": "vladsaveliev"
    }
  ],
  "successCodes": [],
  "label": "TargQC",
  "sbg:id": "vladsaveliev/testproject/targqc/24",
  "sbg:revision": 24,
  "arguments": [
    {
      "separate": false,
      "valueFrom": {
        "class": "Expression",
        "script": "{\n  return \"&& python html_to_b64html.py ./targqc_results/summary.html\"\n}",
        "engine": "#cwl-js-engine"
      },
      "position": 999
    },
    {
      "separate": true,
      "valueFrom": "./targqc_results",
      "prefix": "-o",
      "position": -1
    }
  ],
  "sbg:modifiedOn": 1472567627,
  "inputs": [
    {
      "type": [
        {
          "type": "array",
          "name": "bams",
          "items": "File"
        }
      ],
      "label": "BAM files",
      "inputBinding": {
        "separate": true,
        "secondaryFiles": [
          ".bai"
        ],
        "sbg:cmdInclude": true,
        "position": 0,
        "itemSeparator": " "
      },
      "description": "BAM files for evaluation",
      "sbg:stageInput": "link",
      "sbg:fileTypes": "BAM,FQ,FASTQ,FQ.GZ,FASTQ.GZ",
      "id": "#bams"
    },
    {
      "type": [
        "null",
        "File"
      ],
      "label": "Target BED file",
      "inputBinding": {
        "separate": true,
        "sbg:cmdInclude": true,
        "loadContents": true,
        "prefix": "--bed",
        "position": 2
      },
      "description": "Target BED file",
      "sbg:stageInput": "link",
      "sbg:fileTypes": "BED",
      "id": "#bed"
    },
    {
      "type": [
        "null",
        "int"
      ],
      "label": "Threads",
      "inputBinding": {
        "separate": true,
        "sbg:cmdInclude": true,
        "prefix": "-t"
      },
      "description": "Threads number",
      "sbg:toolDefaultValue": "1",
      "id": "#threads"
    }
  ],
  "sbg:validationErrors": [],
  "id": "https://api.sbgenomics.com/v2/apps/vladsaveliev/testproject/targqc/24/raw/",
  "stdin": "",
  "sbg:job": {
    "inputs": {
      "threads": null,
      "bed": {
        "size": 0,
        "secondaryFiles": [],
        "class": "File",
        "path": "/path/to/target.bed"
      },
      "bams": [
        {
          "size": 0,
          "secondaryFiles": [],
          "class": "File",
          "path": "/path/to/sample-1.bam"
        },
        {
          "size": 0,
          "secondaryFiles": [],
          "class": "File",
          "path": "/path/to/sample-2.bam"
        }
      ]
    },
    "allocatedResources": {
      "cpu": 1,
      "mem": 1000
    }
  },
  "sbg:cmdPreview": "targqc -o ./targqc_results  /path/to/sample-1.bam /path/to/sample-2.bam && python html_to_b64html.py ./targqc_results/summary.html",
  "sbg:sbgMaintained": false,
  "requirements": [
    {
      "class": "ExpressionEngineRequirement",
      "id": "#cwl-js-engine",
      "requirements": [
        {
          "class": "DockerRequirement",
          "dockerPull": "rabix/js-engine"
        }
      ]
    },
    {
      "class": "CreateFileRequirement",
      "fileDef": [
        {
          "filename": "html_to_b64html.py",
          "fileContent": "import os\nimport sys\nimport base64\nimport mimetypes\nfrom bs4 import BeautifulSoup\n\n\ndef execute(html_fpath):\n    b64_html_fpath = os.path.splitext(html_fpath)[0] + '.b64html'\n    if os.path.isfile(html_fpath):\n        with open(b64_html_fpath, 'wt') as fp:\n            fp.write(html_to_dataurl(html_fpath))\n\n\ndef dataurl(data, mime=None):\n    isfile = os.path.isfile(data)\n    if not isfile and not mime:\n        raise Exception('Mimetype must be provided when encoding data is not a valid file path.')\n    if not mime:\n        mimetypes.init()\n        mime, enc = mimetypes.guess_type(os.path.join('file://', data))\n        if mime is None:\n            raise Exception('rfc2397: failed to determine file type')\n    if isfile:\n        with open(data, 'r') as fp:\n            data = fp.read()\n    return 'data:%s;base64,%s' % (mime, base64.b64encode(data))\n\n\ndef compact_html(html_file):\n    with open(html_file) as f:\n        html = f.read()\n    base_dir = os.path.split(html_file)[0]\n    soup = BeautifulSoup(html)\n    for img in soup.findAll('img'):\n        durl_img = dataurl(os.path.join(base_dir, img['src']))\n        img['src'] = durl_img\n\n    js = \"javascript: void(0); document.getElementById('%s').scrollIntoView(true);\"\n    for anchor in soup.findAll('a'):\n        if 'href' in anchor:\n            continue\n        if anchor['href'].startswith('#'):\n            anchor['href'] = js % anchor['href'][1:]\n        else:\n            del anchor['href']\n    return soup.prettify()\n\n\ndef html_to_dataurl(html_file):\n    return dataurl(compact_html(html_file), mime='text/html')\n\n\nif __name__ == \"__main__\":\n    execute(sys.argv[1] if len(sys.argv) > 1 else 'targqc_results/summary.html')"
        }
      ]
    }
  ]
}