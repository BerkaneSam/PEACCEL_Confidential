import openpyxl
import argparse


def get_arguments():
    parser = argparse.ArgumentParser(description="Splits xlsx sheets into different files")
    parser.add_argument('data', metavar='data', type=str, help="path to the xlsx file")
    return parser.parse_args()


def dataset_making(data):
    """
    Retrieving data from a given xlsx sheet
    :param data: xlsx sheet
    :return: list of data
    """
    datas = []
    for row in data.values:
        datas.append(row)
    return datas


def get_data(path):
    """
    Retrieves xlsx data contained in sheets into a list, separates the sheets
    :param path: path to xlsx files
    :return:
    """
    datasets = []
    iterated_sheet = 1
    workbook = openpyxl.load_workbook(path)
    sheets = workbook.sheetnames
    sheets_num = len(workbook.worksheets)
    for sname in sheets:
        datasets.append([sname, dataset_making(workbook[sname])])
        print(f"  sheet {iterated_sheet}/{sheets_num} treated")
        iterated_sheet += 1
    return datasets


def data_treatment_making(datasets):
    """
    Takes all data contained in different xlsx sheets and makes new csv files(datasets) each corresponding to a given
    sheet
    :param datasets: list of xlsx sheets
    :return: makes multiple csv files
    """
    data_size = len(datasets)
    count = 1
    print(data_size)
    for data in datasets:
        with open(f"{data[0]}.csv", 'w')as filout:
            print(data[0])
            for i in range(len(data[1])):
                line_to_write = ""
                for j in range(len(data[1][i])):
                    if i == 0 and j == 0:
                        continue
                    elif j == 0:
                        line_to_write += str(data[1][i][j])
                    elif data[1][i][j] is None:
                        line_to_write = line_to_write + ","
                    else:
                        line_to_write = line_to_write + "," + str(data[1][i][j])
                filout.write(f"{line_to_write}\n")
        print(f"  iteration {count}/{data_size} treated")
        count += 1


if __name__ == '__main__':
    print("program launched")
    args = get_arguments()
    print("arguments retrieved")
    print("data pretreatment ongoing...")
    pretreated_data = get_data(args.data)
    print("data pretreatment done")
    print("data treatment ongoing...")
    data_treatment_making(pretreated_data)
    print("program done")
