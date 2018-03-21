import pyfits
import os
import re
from collections import OrderedDict


def JD_form_fits_files(dir_path):
    fit_files = [file for file in os.scandir(dir_path) if file.name.endswith('.fits')]

    m_date = OrderedDict()
    m_frame = list()
    for file in sorted(fit_files, key=lambda x: x.name):
        file_frame, file_filt = re.findall(r'.*-(\d+)_(\w+)\.fits$', file.path)[0]
        with pyfits.open(file.path) as f:
            JD = str(f[0].header['JD'])
        if file_frame not in m_frame:
            m_frame.append(file_frame)
        if file_filt in m_date:
            m_date[file_filt][file_frame] = JD
        else:
            m_date[file_filt] = {file_frame: JD}

    save_path = os.path.join(dir_path, 'frameJD.txt')
    with open(save_path, 'w') as f:
        f.write('Frame,{}\n'.format(','.join(m_date.keys())))
        for frame in sorted(m_frame):
            temp = ['{}'.format(frame)]
            for filt in m_date.keys():
                if frame in m_date[filt]:
                    temp.append(m_date[filt][frame])
                else:
                    temp.append('-')
            f.write('{}\n'.format(','.join(temp)))
    return save_path
