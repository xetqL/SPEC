import argparse
import numpy

def flatten(l):
    return [item for sublist in l for item in sublist] #just flatten once


if __name__ == '__main__':
    print("===========================")
    print("  Output File Parser v1.0  ")
    print("===========================")

    parser = argparse.ArgumentParser(description='Process the outputfile.')
    parser.add_argument('--std',  type=str, help='the outputfile std  to process.')
    parser.add_argument('--ulba', type=str, help='the outputfile ulba to process.')
    args = parser.parse_args()

    std_of  = open(args.std, 'r')
    content = std_of.read()
    lines   = content.split('\n')
    total_time_std = float([l.split(':')[-1] for l in lines if 'total_time' in l][0])
    print("Total time for the standard approach is:", total_time_std)
    avg_step_time_std = numpy.mean(flatten([map(float, l.split(':')[-1].split(', ')) for l in lines if '"step times"' in l]))
    avg_comp_time_std = numpy.mean(flatten([map(float, l.split(':')[-1].split(', ')) for l in lines if '"comp times"' in l]))

    ulba_of  = open(args.ulba, 'r')
    content = ulba_of.read()
    lines   = content.split('\n')
    total_time_ulba = float([l.split(':')[-1] for l in lines if 'total_time' in l][0])
    print("Total time for the ULBA approach is: \t", total_time_ulba)

    avg_step_time_ulba = numpy.mean(flatten([map(float, l.split(':')[-1].split(', ')) for l in lines if '"step times"' in l]))
    avg_comp_time_ulba = numpy.mean(flatten([map(float, l.split(':')[-1].split(', ')) for l in lines if '"comp times"' in l]))

    print("\nAverage step time (std)  {0:.3f} s".format(avg_step_time_std))
    print(  "Average step time (ulba) {0:.3f} s".format(avg_step_time_ulba))
    print("\nAverage comp time (std)  {0:.3f} s".format(avg_comp_time_std))
    print(  "Average comp time (ulba) {0:.3f} s".format(avg_comp_time_ulba))

    print("\nPerformance improvement: {0:.3f} %".format(100.*(total_time_std-total_time_ulba)/total_time_std))

