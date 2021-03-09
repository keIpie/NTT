#!/usr/bin/env sage

import argparse, os, sys, time

from ntt import *

def main():
    parser = argparse.ArgumentParser(description='PRF representations.')
    parser.add_argument("-ntt", help="test NTT implementation", action="store_true")
    parser.add_argument("-n", type=int, help="cyclotomic polynomial power", action="store", default=1024)
    parser.add_argument("-q", type=int, help="modulus", action="store", default=12289)
    parser.add_argument("-verbose", help="activate printing additional information", action="store_true")
    args = parser.parse_args()

    ##########################################################################
    #      TESTING NTT IMPLEMeNTATION
    ##########################################################################
    if args.ntt:
        test_ntt()


if __name__ == "__main__":
    main()
