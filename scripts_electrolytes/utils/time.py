from datetime import datetime, timedelta

def when_is_now():
    ''' Simply formatting the current date and time'''
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S")


def increase_jobtime(timestr, factor):
    ''' Increase the Slurm time request and format it correctly'''
    split = timestr.split(':')

    if len(split) == 3:
        # Only HH:MM:SS or H:MM:SS
        hh, mm, ss = split
        dd = 0

    elif len(split) == 4:
        # D:HH:MM:SS or DD:HH:MM:SS
        dd, hh, mm, ss = split


    seconds = 24*3600*dd + 3600*hh + 60*mm + ss
    delta = datetime.timedelta(seconds=int(factor*seconds))

    try:
        d = str(delta).split(' days')[0]
        if d>=7:
            raise ValueError('Trying to launch a job of more than 7 days. Aborting,')

        h,m,s = str(delta).split(', ')[1].split(':')
        h = '{0:02d}'.format(int(h))
        time = '{}:{}:{}:{}'.format(d, h,m,s)
    except:
        h,m,s = str(delta).split(':')
        time = '{}:{}:{}'.format(h,m,s)

    return time
