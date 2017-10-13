#! /usr/bin/env python
"""
Script to combine two makesourcedb sky models
"""
import argparse
from argparse import RawTextHelpFormatter
import lsmtool
import sys
import os


def main(model_list, skymodel):
    """
    Combines makesourcedb sky models

    Parameters
    ----------
    model_list : str
        Filenames of the input makesourcedb sky model 1
    skymodel : str
        Filename of the output makesourcedb sky model

    """
    model_list = model_list.strip('[]').split(',')
    model_list  = [f.strip() for f in model_list]

   # First find a model with sources
    for sm in model_list[:]:
        model_list.remove(sm)
        try:
            s1 = lsmtool.load(sm)
            break
        except:
            pass

    # Now try to load the rest of the sky models and combine with first one
    for sm in model_list:
        try:
            s2 = lsmtool.load(sm)

            # Combine sky models, keeping all sources
            s1.ungroup()
            s2.ungroup()
            s1.concatenate(s2, keep='all')
        except:
            pass

    s1.write(skymodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Combine two makesourcedb sky models.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('model1', help='name of input makesourcedb sky model 1')
    parser.add_argument('model2', help='name of input makesourcedb sky model 2')
    parser.add_argument('skymodel', help='name of the output makesourcedb sky model')
    args = parser.parse_args()

    main(args.model1, args.model2, args.skymodel)
