import argparse
import io

def main(output_file_name, names_list, numeric_list, image_list, log_list):

    output_file = io.open(output_file_name, 'w', encoding='utf8')
    
    names_array = []
    if names_list is not None:
        with open(names_list, "r") as names_f: 
            names_array = names_f.read().splitlines()
    
    
    if numeric_list is not None:
        with open(numeric_list, "r") as numeric_f: 
            numeric_array = numeric_f.read().splitlines()
    
    
    if image_list is not None:
        with open(image_list, "r") as image_f: 
            image_array = image_f.read().splitlines()

    
    if log_list is not None:
        with open(log_list, "r") as log_f: 
            log_array = log_f.read().splitlines()


    csv_values = numeric_array + image_array + log_array
    for idx in range(len(names_array)):
        output_file.write(names_array[idx] + ", " + csv_values[idx] + "\n")








if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)
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