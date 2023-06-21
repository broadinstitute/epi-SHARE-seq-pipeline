import argparse
import io

def main(output_file_name, names_list, numeric_list, image_list, log_list):

    output_file = io.open(output_file_name, 'w', encoding='utf8')
    
    if names_list is not None:
        with open(names_list, "r") as names_f: 
            names_array = names_f.readline()

    if numeric_list is not None:
        with open(numeric_list, "r") as names_f: 
            numeric_array = names_f.readline()

    if image_list is not None:
        with open(names_list, "r") as names_f: 
            names_array = names_f.readline()

    if log_list is not None:
        with open(names_list, "r") as names_f: 
            log_array = names_f.readline()

    csv_values = names_list.append(numeric_list.append(image_list.append(log_list)))
    for idx in range(len(names_list)):
        output_file.write(names_list[idx] + " " + csv_values[idx] + ",")








if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('output_file_name',
                       help='html file to write to')
    group.add_argument('names_list', 
                       help='file containing list of names of relevant data from a pipeline run')
    group.add_argument('numeric_list',
                       help='file containing list of numbers produced by the pipeline')
    group.add_argument('image_list',
                       help='file containing list of png strings?') 
    
    group.add_argument('log_list',
                       help='file containing the log information for a run')
    args = parser.parse_args()
    main(args.output_file_name, args.names_list, args.numeric_list, args.image_list, args.log_list)