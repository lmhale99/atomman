from pathlib import Path
import os
import glob

notebookdir = Path('tutorial')
docdir = Path('source', 'tutorial')

# Delete old Notebooks from the source folder
for notebook in docdir.glob('*.ipynb'):
    notebook.unlink()

pagenames = []

# Copy new Notebooks from tutorial to source
for notebook in notebookdir.glob('*.ipynb'):
    pagenames.append(notebook.stem)
    # Read in Notebook content
    with open(notebook) as f:
        content = f.read()
    
    # Change all links from .ipynb to .html
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
    oldname = notebook.name
    newname = oldname.replace(' ', '_')
    newpath = Path(docdir, newname)
    with open(newpath, 'w') as f:
        f.write(content)

# Define sortkey for page numbers
def sortkey(fname):
    """
    Sorts filenames based on numbers at the beginning of the file with up to
    two subsets, i.e. #., #.#., or #.#.#.
    """
    numbers = fname.split(' ')[0].split('.')
    val = int(numbers[0]) * 10000
    try:
        val += int(numbers[1]) * 100
    except:
        pass
    else:
        try:
            val += int(numbers[2])
        except:
            pass
    return val

# Create index.rst header
indexlines = [
    "Tutorials",
    "=========",
    '',
    '.. toctree::',
    '   :maxdepth: 2',
    '   ']

# Build lines for each Notebook
for name in sorted(pagenames, key=sortkey):
    indexlines.append(f'   {name} <{name.replace(" ", "_")}>')

with open(Path(docdir, 'index.rst'), 'w', encoding='UTF-8') as f:
    f.write('\n'.join(indexlines))