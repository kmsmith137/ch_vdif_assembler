#!/usr/bin/env python

import os
import re
import sys


def index_thumbnails(dirname):
    imgdir = os.path.join(dirname, 'img')

    fname_hash = { }   # (imtype,pol,frame_ix) -> filename
    all_indices = set()

    for fname in os.listdir(imgdir):
        for imtype in [ 'full', 'thumbnail' ]:
            for pol in [ 0, 1 ]:
                prefix = '%s_pol%d_frame' % (imtype, pol)
                if not fname.startswith(prefix):
                    continue
            
                m = re.match(r'(\d+)\.png', fname[len(prefix):])
                if not m:
                    continue

                frame_ix = int(m.group(1))
                fname_hash[imtype,pol,frame_ix] = fname
                all_indices.add(frame_ix)

    
    html_filename = os.path.join(dirname, 'index.html')
    html = open(html_filename, 'w')

    print >>html, '<html><head><title>%s</title></head>' % os.path.basename(dirname)
    print >>html, '<body><table cellspacing="10">'

    for frame_ix in sorted(all_indices):
        print >>html, '<tr> <td> Frame %d' % frame_ix

        for pol in [0,1]:
            full = fname_hash.get(('full',pol,frame_ix), None)
            thumb = fname_hash.get(('thumbnail',pol,frame_ix), None)

            if thumb and full:
                print >>html, '<td> <a href="img/%s"> <img src="img/%s"> </a>' % (full, thumb)
            elif thumb and not full:
                print >>html, '<td> <img src="img/%s">' % thumb
            elif not thumb and full:
                print >>html, '<td> <a href="img/%s"> pol %d </a>' % (full, pol)
            else:
                print >>html, '<td>'

    print >>html, '</table></body></html>'
    print >>sys.stderr, 'wrote', html_filename


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print >>sys.stderr, 'usage: index-thumbnails.py <waterfall_outdir>'
        sys.exit(2)

    index_thumbnails(sys.argv[1])
