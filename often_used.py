import csv


def open_csv(csv_file_path, delimiter=',', new_line=False):
    '''
    Open XMatch_NOMAD_USNO_path
    :param csv_file_path: Path to XMatch NOMAD1 and USNO-B1 catalogues
    :param delimiter: delimiter in csv-file
    :param new_line: Skip (True) the first line or not (False)
    :return: data of stars in XMatch_NOMAD_USNO_path file
    '''

    with open(csv_file_path, 'r') as f:
        header = f.readline().strip().split(delimiter)

    with open(csv_file_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter=delimiter)
        if new_line:
            next(csv_file)
        stars_from_file = [row for row in csv_file]

    return header, stars_from_file


def write_to_csv(save_path, data, header, delimiter=','):
    with open(save_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(data)
    print('File {} has created!'.format(save_path))