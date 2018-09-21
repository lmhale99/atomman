from __future__ import print_function, unicode_literals
import os
import glob

for oldpath in glob.iglob(os.path.join('tutorial', '*.ipynb')):
    
    # Read in Notebook content
    with open(oldpath) as f:
        content = f.read()
    
    # Change all links to Notebooks to html links
    for i in range(100):
        try:
            endindex = content.index('.ipynb)')+6
        except:
            break
        else:
            startindex = content[:endindex].rindex('](') + 2
            oldlink = content[startindex:endindex]
            newlink = oldlink.replace('.ipynb', '.html').replace(' ', '_')
            content = content.replace(oldlink, newlink)
    
    # Save to source directory
    oldname = os.path.basename(oldpath)
    newname = oldname.replace(' ', '_')
    newpath = os.path.join('source', 'tutorial', newname)
    with open(newpath, 'w') as f:
        f.write(content)