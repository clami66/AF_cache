import json
import argparse
from string import ascii_uppercase, ascii_lowercase
ascii_upperlower = ascii_uppercase + ascii_lowercase
from alphafold3.data.msa import Msa

def main(args):
    parsed_jsons = [json.load(open(j)) for j in args.jsons]

    main_json = parsed_jsons[0]

    for i, extra_json in enumerate(parsed_jsons[1:]):
        chain = ascii_upperlower[i+1]
        sequences = extra_json["sequences"]
        sequences[0]["protein"]["id"] = chain
        main_json["sequences"].append(sequences[0])

    msapath = main_json["sequences"][0]["protein"]["unpairedMsaPath"]
    msaquery = main_json["sequences"][0]["protein"]["sequence"]
    with open(msapath) as m:
        #msa = Msa.from_a3m(msaquery, "protein", m.readlines())
        a3m = "".join(m.readlines())

    merged_msa = Msa.from_multiple_a3ms([a3m, a3m], "protein", False)
    print(merged_msa.to_a3m())
    with open(args.out_json, "w") as out:
      json.dump(main_json, out, indent=4, sort_keys=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merges multiple inputs to AF3 (json files) into a single, multichain input. Chains are named according to the order of json files submitted by the user")
    parser.add_argument("--jsons", nargs = "+", help = "Path to multiple AF3 inputs (json files)")
    parser.add_argument("--out_json", help = "Path to output json file")
    args = parser.parse_args()

    main(args)
