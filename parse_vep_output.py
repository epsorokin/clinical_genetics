
import argparse
import re
import sys

def main(args):

    f = open(args.vcf)
    extra_fields = None
    header = None
    for line in f:
        line=line.strip()
        if line.startswith('#'):
            line=line.lstrip('#')
            if line.find('IMPACT=')>0:
                extra_fields=line.split('\t')[-1].split(';')
        continue

    fields=line.split('\t')
    
    try: 
        print '%s\t%s' % (line.split('\t')[0], extra_fields[args.field] ) 
    except KeyError:
        print '%s\t%s' % (line.split('\t')[0], 'NA')

    if extra_fields is None:
        print >> sys.stderr, "Error: No annotation fields to parse."
        sys.exit(1)
    f.close()

if __name__== '__main__':

    myparser = argparse.ArgumentParser()
    myparser.add_argument ('--vcf',required=True)
    myparser.add_argument ('--field',help="Select annotation field",required=True)
    args = myparser.parse_args()
    main(args)
