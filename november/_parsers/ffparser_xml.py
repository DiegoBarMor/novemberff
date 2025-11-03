from enum import Enum, auto

import november as nov

# ////////////////////////////////////////////////////////////////////////////////
class State(Enum):
    START = auto()
    TAG_NAME = auto()
    TAG_CONTENT = auto()
    ATTRIBUTE_KEY = auto()
    ATTRIBUTE_VAL = auto()
    END = auto()


# ////////////////////////////////////////////////////////////////////////////////
class XMLNode:
    def __init__(self, tag_name, attributes, tag_content):
        self.tag_name: str = tag_name
        self.attributes: dict[str, str] = attributes
        self.tag_content: str = tag_content
        self.children: list["XMLNode"] = []

    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"XMLNode(name='{self.tag_name}', nattrs={len(self.attributes)}, content='{self.tag_content}', nchildren={len(self.children)})"

    # --------------------------------------------------------------------------
    def add_child(self, node):
        self.children.append(node)

    # --------------------------------------------------------------------------
    def get_child_by_name(self, name: str) -> "XMLNode | None":
        for child in self.children:
            if child.tag_name == name:
                return child

    # --------------------------------------------------------------------------
    def has_attr(self, key: str) -> bool:
        return key in self.attributes.keys()

    # --------------------------------------------------------------------------
    def get_attr(self, key: str) -> str | None:
        if self.has_attr(key):
            return self.attributes[key]


# ////////////////////////////////////////////////////////////////////////////////
class FFParserXML(nov.ForceFieldParser):
    def __init__(self, path_xml):
        self.root = XMLNode('', {}, '')
        self._path_xml = path_xml
        with open(path_xml, 'r') as file:
            self._xml_content = file.read()


    # --------------------------------------------------------------------------
    def parse(self):
        node = self.root
        parent_chain = []
        for line in self._xml_content.splitlines():
            tag_name, attributes, tag_content, ended = self.parse_xml_line(line)

            ### CASE 0: end-tag only
            if not tag_name:
                node = parent_chain.pop()
                continue

            child = XMLNode(tag_name, attributes, tag_content)

            ### CASE 1: empty-element-tag or (start-tag + content + end-tag) in a single line
            node.add_child(child)
            if ended: continue

            ### CASE 2: start-tag only
            parent_chain.append(node)
            node = child


    # --------------------------------------------------------------------------
    @classmethod
    def parse_xml_line(cls, line: str):
        state = State.START
        tag_name = ""
        tag_content = ""
        attributes = {}
        current_key = ""

        line = line.strip()
        for i,c in enumerate(line):
            if state == State.START:
                if line[i+1] == "/":
                    state = State.END
                elif c == "<":
                    state = State.TAG_NAME
                continue

            if state == State.TAG_NAME:
                if c == ">":
                    state = State.TAG_CONTENT
                elif c == " ":
                    state = State.ATTRIBUTE_KEY
                else:
                    tag_name += c
                continue

            if state == State.TAG_CONTENT:
                if c == "<":
                    state = State.END
                else:
                    tag_content += c
                continue

            if state == State.ATTRIBUTE_KEY:
                if c == "=":
                    attributes[current_key] = ""
                elif c == "\"":
                    state = State.ATTRIBUTE_VAL
                elif c == ">":
                    state = State.TAG_CONTENT
                elif c == "/":
                    state = State.END
                elif c != " ":
                    current_key += c
                continue

            if state == State.ATTRIBUTE_VAL:
                if c == "\"":
                    current_key = ""
                    state = State.ATTRIBUTE_KEY
                elif c != " ":
                    attributes[current_key] += c
                continue

        return tag_name, attributes, tag_content, state == State.END


# ////////////////////////////////////////////////////////////////////////////////
