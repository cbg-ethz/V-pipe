import re
import os
import argparse

def extract_code_chunks(markdown_text):
    """
    Extracts code chunks from markdown text, specifically those starting with 
    ```bash and excluding chunks with '## not run' on the second line.

    Args:
        markdown_text: The markdown text to process.

    Returns:
        A list of code blocks.
    """
    pattern = r"```bash\n(?!## not run\n)(.*?)\n```"
    return re.findall(pattern, markdown_text, re.DOTALL)

def process_markdown_file(filename, output_dir):
    """
    Reads a markdown file, extracts bash code chunks, and saves them to a 
    new file with the .sh extension in the specified output directory.

    Args:
        filename: The path to the markdown file.
        output_dir: The directory where the output .sh files will be saved. 
                   Defaults to the current directory.
    """
    try:
        with open(filename, 'r') as f:
            markdown_content = f.read()

        code_blocks = extract_code_chunks(markdown_content)

        if code_blocks:
            output_filename = os.path.join(output_dir, os.path.basename(os.path.splitext(filename)[0]) + '.sh')
            with open(output_filename, 'w') as f:
                for code_block in code_blocks:
                    f.write(code_block + '\n')
            print(f"Bash code extracted to: {output_filename}")
        else:
            print(f"No bash code blocks found in: {filename}")

    except FileNotFoundError:
        print(f"File not found: {filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract bash code from markdown files.")
    parser.add_argument("markdown_files", nargs="+", help="Paths to the markdown files.")
    parser.add_argument("-o", "--output_dir", default=os.getcwd(), 
                        help="Directory to save the output .sh files. Defaults to the current directory.")
    args = parser.parse_args()

    for markdown_file in args.markdown_files:
        process_markdown_file(markdown_file, args.output_dir)