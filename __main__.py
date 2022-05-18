import argparse

from utils import classify_mutation

parser = argparse.ArgumentParser()
parser.add_argument('-normal')
parser.add_argument('-mutant')
args = parser.parse_args()

def main(normal, mutant):
    return classify_mutation(normal, mutant)

if __name__ == "__main__":
    main(args.normal, args.mutant) 
