#pragma once
#include <cstddef>
#include <vector>
#include <unordered_set>
#include <queue>
#include "boundingbox.h"
#include <fstream>

/*!
 * \brief Шаблонный класс RStarTree
 * \param[in] LeafType - тип, хранимый в листе дерева
 * \param[in] dimensions - количество измерений огранич. области
 * \param[in] min_child_items - минимальное количество детей узла
 * \param[in] max_child_items - максимальное количество детей узла
 *  */
template<
    typename LeafType,
    std::size_t dimensions,
    std::size_t min_child_items,
    std::size_t max_child_items
>
class RStarTree {
    using BoundingBox = RStarBoundingBox<dimensions>;
private:

	
    struct TreePart {
        BoundingBox box;  
    };

	
    struct Node : public TreePart {
        std::vector<TreePart*> items;
        int hasleaves{ false };
    };

	
    struct Leaf : public TreePart {
        LeafType value;
    };

    struct SplitParameters {
        int index{ -1 };
        int axis{ -1 };
        axis_type type{ axis_type::lower };
    };

public:
    //using BoundingBox = RStarBoundingBox<dimensions>;

    /*!
	* \brief Лист, экземпляры которого находятся в массиве,
	*  возвращаемом в методе поиска, с константной ограничивающей областью,
	*  обращение к полям посредством методов
	* */
    class LeafWithConstBox { 
    public:                  
        LeafWithConstBox(Leaf* leaf_) :leaf(leaf_) {}
        /*!
		* \brief Возвращает константную ссылку на огранич. область
		* */
        const BoundingBox& get_box() const {
            return leaf->box;
        }
        const BoundingBox& get_box() {
            return leaf->box;
        }
        const LeafType& get_value() const {
            return leaf->value;
        }
        /*!
		* \brief Возвращает ссылку на значение листа
		* */
        LeafType& get_value() {
            return leaf->value;
        }
        bool operator<(const LeafWithConstBox& rhs) const {
            if (get_value() == rhs.get_value()) {
                return get_box() < rhs.get_box();
            }
            return get_value() < rhs.get_value();
        }
        bool operator==(const LeafWithConstBox& rhs) const {
            return get_value() == rhs.get_value() && get_box() == rhs.get_box();
        }
        bool operator!=(const LeafWithConstBox& rhs) const {
            return !operator==(rhs);
        }
    private:
        Leaf* leaf{ nullptr };
    };


public:

    //! \brief Конструктор RStarTree, проверяет шаблонные параметры на корректность
    RStarTree() {
        if (dimensions <= 0 || max_child_items < min_child_items) {
            throw std::invalid_argument("");
        }
    }
    ~RStarTree() {
        delete_tree(tree_root);
    }

    /*!
 	* \brief Метод вставки листа в дерево
 	* \param[in] leaf - значение листа
 	* \param[in] box - ограничивающая область
 	*  */
    void insert(const LeafType leaf, const BoundingBox& box) {
        size_++;
        Leaf* new_leaf = new Leaf;
        new_leaf->value = leaf;
        new_leaf->box = box;
        if (!tree_root) {
            tree_root = new Node();
            tree_root->hasleaves = true;
            tree_root->items.reserve(min_child_items);
            tree_root->items.push_back(new_leaf);
        }
        else {
            choose_leaf_and_insert(new_leaf, tree_root);
        }
        used_deeps.clear();
    }

    /*!
 	* \brief Метод поиска объектов в области, возвращается массив листьев
 	* \param[in] box - область поиска
 	*  */
    std::vector<LeafWithConstBox> find_objects_in_area(const BoundingBox& box) {
        std::vector<LeafWithConstBox>leafs;
        find_leaf(box, leafs, tree_root);
        return leafs;
    }

    /*!
 	* \brief Метод удаления объектов в области
 	* \param[in] box - область поиска
 	*  */
    void delete_objects_in_area(const BoundingBox& box) {
        delete_leafs(box, tree_root);

    }

    /*!
 	* \brief Метод записи дерева в бинарном формате
 	* \param[in] file - поток, в который ведется запись дерева
 	*  */
    void write_in_binary_file(std::fstream& file) {
        int dimensions_ = dimensions;
        int min_child_items_ = min_child_items;
        int max_child_items_ = max_child_items;
        file.write(reinterpret_cast<const char*>(&dimensions_), sizeof(dimensions_));
        file.write(reinterpret_cast<const char*> (&min_child_items_), sizeof(min_child_items_));
        file.write(reinterpret_cast<const char*>(&max_child_items_), sizeof(max_child_items_));
        std::queue<Node*> nodes;//Ñ ïîìîùüþ BFS ïîóðîâíåâî âûâîäÿòñÿ âåðøèíû â ôàéë
        nodes.push(tree_root);
        while (!nodes.front()->hasleaves) {//Âûâåäóòñÿ âåðøèíû, äåòè êîòîðûõ íå ëèñòüÿ
            write_node(nodes.front(), file);
            for (size_t i = 0; i < nodes.front()->items.size(); i++) {
                nodes.push(static_cast<Node*>(nodes.front()->items[i]));
            }
            nodes.pop();
        }
        std::vector<Leaf*>leafs;
        while (!nodes.empty()) {
            for (long long i = 0; i < nodes.front()->items.size(); i++) {
                leafs.push_back(static_cast<Leaf*>(nodes.front()->items[i]));// Ñîõðàíÿþòñÿ â ìàññèâ ëèñòüÿ
            }                                                                // ïîñëå âûâîäèòñÿ ñàìà âåðøèíà
            write_node(nodes.front(), file);
            nodes.pop();
        }
        for (auto& leaff : leafs) {
            write_leaf(leaff, file);// Âûâîäÿòñÿ ëèñòüÿ
        }
    }

    /*!
 	* \brief Метод чтения дерева в бинарном формате
 	* \param[in] file - поток, из которого ведется чтение дерева
 	*  */
    void read_from_binary_file(std::fstream& file) {
        int dimensions_ = dimensions;
        int min_child_items_ = min_child_items;
        int max_child_items_ = max_child_items;
        file.read(reinterpret_cast<char*>(&dimensions_), sizeof(dimensions_));
        file.read(reinterpret_cast<char*> (&min_child_items_), sizeof(min_child_items_));
        file.read(reinterpret_cast<char*>(&max_child_items_), sizeof(max_child_items_));
        if (dimensions_ != dimensions || min_child_items != min_child_items_ || max_child_items != max_child_items_)
            throw std::invalid_argument("");
        if (tree_root) {
            delete_tree(tree_root);//óäàëÿåì äåðåâî, åñëè îíî áûëî
        }
        tree_root = read_node(file);
        std::queue<Node*> nodes;// Íàîáîðîò ñ ïîìîùüþ BFS ñ÷èòûâàþòñÿ âåðøèíû èç ôàéëà
        nodes.push(tree_root);
        while (!nodes.empty()) {
            for (size_t i = 0; i < nodes.front()->items.size(); i++) {
                Node* new_node = read_node(file);
                nodes.front()->items[i] = new_node;
                nodes.push(new_node);
            }
            nodes.pop();
            if (nodes.front()->hasleaves) break;
        }
        while (!nodes.empty()) {
            for (size_t i = 0; i < nodes.front()->items.size(); i++) {
                Leaf* new_leaf = read_leaf(file);
                nodes.front()->items[i] = new_leaf;
            }
            nodes.pop();
        }
    }
private:
    void write_node(Node* node, std::fstream& file) {
        size_t _size = node->items.size();
        file.write(reinterpret_cast<const char*>(&_size),
            sizeof(_size));
        file.write(reinterpret_cast<const char*>(&(node->hasleaves)), sizeof(node->hasleaves));
        for (size_t axis = 0; axis < dimensions; axis++) {
            file.write(reinterpret_cast<const char*>(&(node->box.max_edges[axis])),
                sizeof(node->box.max_edges[axis]));
            file.write(reinterpret_cast<const char*>(&(node->box.min_edges[axis])),
                sizeof(node->box.min_edges[axis]));
        }
    }

    void write_leaf(Leaf* leaf, std::fstream& file) {
        for (size_t axis = 0; axis < dimensions; axis++) {
            file.write(reinterpret_cast<char*>(&(leaf->box.max_edges[axis])),
                sizeof(leaf->box.max_edges[axis]));
            file.write(reinterpret_cast<char*>(&(leaf->box.min_edges[axis])),
                sizeof(leaf->box.min_edges[axis]));
        }
        file.write(reinterpret_cast<char*>(&(leaf->value)),
            sizeof(LeafType));
        size_--;
    }

    Node* read_node(std::fstream& file) {
        size_t _size;
        file.read(reinterpret_cast<char*>(&_size), sizeof(_size));
        Node* new_node = new Node;
        new_node->items.resize(_size);
        file.read(reinterpret_cast<char*>(&(new_node->hasleaves)), sizeof(new_node->hasleaves));
        for (size_t axis = 0; axis < dimensions; axis++) {
            file.read(reinterpret_cast<char*>(&(new_node->box.max_edges[axis])),
                sizeof(new_node->box.max_edges[axis]));
            file.read(reinterpret_cast<char*>(&(new_node->box.min_edges[axis])),
                sizeof(new_node->box.min_edges[axis]));
        }
        return new_node;
    }

    Leaf* read_leaf(std::fstream& file) {
        Leaf* new_leaf = new Leaf;
        for (size_t axis = 0; axis < dimensions; axis++) {
            file.read(reinterpret_cast<char*>(&(new_leaf->box.max_edges[axis])),
                sizeof(new_leaf->box.max_edges[axis]));
            file.read(reinterpret_cast<char*>(&(new_leaf->box.min_edges[axis])),
                sizeof(new_leaf->box.min_edges[axis]));
        }
        file.read(reinterpret_cast<char*>(&(new_leaf->value)),
            sizeof(LeafType));
        size_++;
        return new_leaf;
    }

    void find_leaf(const BoundingBox& box, std::vector<LeafWithConstBox>& leafs, Node* node) {
        if (node->hasleaves) {
            for (size_t i = 0; i < node->items.size(); i++) {
                Leaf* temp_leaf = static_cast<Leaf*>(node->items[i]);
                if (box.is_intersected((temp_leaf->box))) {
                    leafs.push_back({ temp_leaf });
                }
            }
        }
        else {
            for (size_t i = 0; i < node->items.size(); i++) {
                Node* temp_node = static_cast<Node*>(node->items[i]);
                if (box.is_intersected((temp_node->box))) {
                    find_leaf(box, leafs, temp_node);
                }
            }
        }
    }

    void delete_leafs(const BoundingBox& box, Node* node) {
        if (node->hasleaves) {//Åñëè äåòè óçëà - ëèñòüÿ, òî ïðîâåðÿþòñÿ âñå äåòè
            for (size_t i = 0; i < node->items.size(); i++) {
                if (box.is_intersected((node->items[i]->box))) {
                    std::swap(node->items[i], node->items.back());// ìåíÿåì èñêîìóþ ñ ïîñëåäíåé
                    delete node->items.back();                    // óäàëÿåì ïîñëåäíþþ
                    node->items.pop_back();
                    i--;
                }
            }
            node->box.reset();
        }
        else {
            for (size_t i = 0; i < node->items.size(); i++) {
                if (box.is_intersected((node->items[i]->box))) {
                    delete_leafs(box, static_cast<Node*>(node->items[i]));// âûçûâàåì ìåòîä îò òåõ âåðøèí,
                }                                                     // îáëàñòè êîòîðûõ ïåðåñåêàþòñÿ ñ èñêîìîé
            }
        }
        for (size_t i = 0; i < node->items.size(); i++) {
            node->box.stretch(node->items[i]->box);
        }

    }

    Node* choose_leaf_and_insert(Leaf* leaf, Node* node, int deep = 0) {
        node->box.stretch(leaf->box);
        if (node->hasleaves) {// Åñëè äåòè óçëà - ëèñòüÿ, òî äîáàâëÿåì ëèñò â ìàññèâ
            node->items.push_back(leaf);
        }
        else {
            Node* new_node = choose_leaf_and_insert(leaf, choose_subtree(node, leaf->box), deep + 1);
            if (!new_node) { // choose_leaf_and_insert âåðíåò óçåë, òî îí âñòàâëÿåòñÿ â ìàññèâ äåòåé
                return nullptr;
            }
            node->items.push_back(new_node);
        }
        if (node->items.size() > max_child_items) {
            return overflow_treatment(node, deep);// Åñëè êîëè÷åñòâî äåòåé áîëüøå max_child_items
        }                                                     // òî óçåë äîëæåí ïîäåëèòüñÿ
        return nullptr;
    }
    Node* choose_node_and_insert(Node* node, Node* parent_node, int required_deep, int deep = 0) {
        parent_node->box.stretch(node->box);
        if (deep == required_deep) {//Àíàëîãè÷åí choose_leaf_and_insert, òîëüêî äëÿ óçëîâ
            parent_node->items.push_back(node);// óçåë óñòàíàâëèâàåòñÿ íà îïðåäåëåííóþ âûñîòó,
        }                                      // ÷òîáû äåðåâî îñòàëîñü èäåàëüíî ñáàëàíñèðîâàííûì
        else {
            Node* new_node = choose_node_and_insert(node, choose_subtree(parent_node, node->box), required_deep, deep + 1);
            if (!new_node) return nullptr;
            parent_node->items.push_back(new_node);
        }
        if (parent_node->items.size() > max_child_items) {
            return overflow_treatment(parent_node, deep);
        }
        return nullptr;
    }
    Node* choose_subtree(Node* node, const BoundingBox& box) {
        std::vector<TreePart*> overlap_preferable_nodes;
        if (static_cast<Node*>(node->items[0])->hasleaves) { //Åñëè äî÷åðíèå óçëû - êîíå÷íûå, 
            int min_overlap_enlargement(std::numeric_limits<int>::max());//èùåòñÿ óçåë ñ íàèìåíüøèì ïåðåêðûòèåì
            int overlap_enlargement(0);
            for (size_t i = 0; i < node->items.size(); i++) {
                TreePart* temp = (node->items[i]);
                overlap_enlargement = box.area() - box.overlap(temp->box);
                if (overlap_enlargement < min_overlap_enlargement) {
                    min_overlap_enlargement = overlap_enlargement;
                    overlap_preferable_nodes.resize(0);
                    overlap_preferable_nodes.push_back(temp);
                }
                else {
                    if (overlap_enlargement == min_overlap_enlargement) {
                        overlap_preferable_nodes.push_back(temp);
                    }
                }
            }
            if (overlap_preferable_nodes.size() == 1) {//Åñëè óçåë îäèí, òî ýòî èñêîìûé óçåë
                return static_cast<Node*>(overlap_preferable_nodes.front());
            }//åñëè íåò, òî èùåòñÿ óçåë, â êîòîðîì áóäåò ìèíèìàëüíîå óâåëè÷åíèå ïëîùàäè
        }
        else {//Åñëè äî÷åðíèå óçëû íå êîíå÷íûå, òî ñîõðàíèì èõ â ìàññèâ
            overlap_preferable_nodes.reserve(node->items.size());
            std::copy(node->items.begin(), node->items.end(), std::back_inserter(overlap_preferable_nodes));
        }
        int min_area_enlargement = std::numeric_limits<int>::max();//è äëÿ êîíå÷íûõ, è äëÿ íå êîíå÷íûõ
        int area_enlargement(0);                                //ïîñëåäóþùèå äåéñòâèÿ îäèíàêîâû
        std::vector<TreePart*>area_preferable_nodes;
        for (size_t i = 0; i < overlap_preferable_nodes.size(); i++) {
            BoundingBox temp(box);
            temp.stretch(overlap_preferable_nodes[i]->box);
            area_enlargement = temp.area() - box.area();
            if (min_area_enlargement > area_enlargement) {
                min_area_enlargement = area_enlargement;
                area_preferable_nodes.resize(0);
                area_preferable_nodes.push_back(overlap_preferable_nodes[i]);
            }
            else {
                if (min_area_enlargement == area_enlargement) {
                    area_preferable_nodes.push_back(overlap_preferable_nodes[i]);
                }
            }
        }
        if (area_preferable_nodes.size() == 1) { //Åñëè óçåë ñ ìèíèìàëüíûì óâåëè÷åíèåì ïëîùàäè îäèí,
            return static_cast<Node*>(area_preferable_nodes.front());//òî âåðíåòñÿ îí
        }
        TreePart* min_area_node{ nullptr };//Èùåòñÿ óçåë ñðåäè îñòàâøèõñÿ ñ ìèíèìàëüíîé ïëîùàäüþ
        int min_area(std::numeric_limits<int>::max());
        for (size_t i = 0; i < area_preferable_nodes.size(); i++) {
            if (min_area > area_preferable_nodes[i]->box.area()) {
                min_area_node = area_preferable_nodes[i];
            }
        }
        return static_cast<Node*>(min_area_node);
    }
    Node* overflow_treatment(Node* node, int deep) {
        if (used_deeps.count(deep) == 0 && tree_root != node) {// Ìåòîä ïîâòîðíîé âñòàâêè âûçûâàåòñÿ
            forced_reinsert(node, deep);                      // òîëüêî îäèí ðàõ íà êàæäîé ãëóáèíå
            return nullptr;                            // è íå äëÿ êîðíÿ
        }
        Node* splitted_node = split(node);// óäàëÿåò ëèøíèõ äåòåé èç node
                                                // è âîçâðàùàåò óçåë ñ íèìè
        if (node == tree_root) {//Åñëè node - êîðåíü äåðåâà,
            Node* temp = new Node;//òî îíà ðàñòåò íà óðîâåíü ââåðõ
            temp->hasleaves = false;
            temp->items.reserve(min_child_items);
            temp->items.push_back(tree_root);
            temp->items.push_back(splitted_node);
            tree_root = temp;
            tree_root->box.reset();
            tree_root->box.stretch(temp->items[0]->box);
            tree_root->box.stretch(temp->items[1]->box);

            return nullptr;
        }
        return splitted_node;//Èíà÷å óçåë âîçâðàùàåòñÿ äëÿ â âñòàâêè â ìàññèâ äåòåé ðîäèòåëÿ

    }
    Node* split(Node* node) {
        SplitParameters params = choose_split_axis_and_index(node);// Âûáèðàåòñÿ íàèáîëåå îïòèìàëüíûå
        std::sort(node->items.begin(), node->items.end(),      // èíäåêñ è îñü
            [&params](auto lhs, auto rhs) {
            return lhs->box.value_of_axis(params.axis, params.type)
                < rhs->box.value_of_axis(params.axis, params.type);
        });
        Node* new_Node = new Node;
        new_Node->items.reserve(max_child_items + 1 - min_child_items - params.index);
        new_Node->hasleaves = node->hasleaves;
        std::copy(node->items.begin() + min_child_items + params.index,
            node->items.end(), std::back_inserter(new_Node->items));
        node->items.erase(node->items.begin() + min_child_items + params.index,
            node->items.end());
        new_Node->box.reset();
        node->box.reset();
        for (auto& w : node->items) {
            node->box.stretch(w->box);
        }
        for (auto& w : new_Node->items) {
            new_Node->box.stretch(w->box);
        }
        return new_Node;
    }
    void forced_reinsert(Node* node, int deep) {//×àñòü äåòåé äàííîé âåðøèíû ïîâòîðíî âñòàâëÿþòñÿ â äåðåâî
        double p = 0.3;//Ïðîöåíò äåòåé, êîòîðûå áóäóò óäàëåíû èç óçëà node
        int number = node->items.size() * p;
        std::sort(node->items.begin(), node->items.end(),
            [&node](auto lhs, auto rhs) {
            return lhs->box.dist_between_centers(node->box) < rhs->box.dist_between_centers(node->box);
        });
        std::vector<TreePart*>forced_reinserted_nodes;
        forced_reinserted_nodes.reserve(number);
        std::copy(node->items.rbegin(), node->items.rbegin() + number, std::back_inserter(forced_reinserted_nodes));
        node->items.erase(node->items.end() - number, node->items.end());
        node->box.reset();
        used_deeps.insert(deep);
        for (TreePart* w : node->items) {
            node->box.stretch(w->box);
        }
        if (node->hasleaves)// åñëè äåòè âåðøèíû - ëèñòüÿ, òî îíè âñòàâëÿþòñÿ ïîâòîðíî ìåòîäîì choose_leaf_and_insert
            for (TreePart* w : forced_reinserted_nodes) {
                choose_leaf_and_insert(static_cast<Leaf*>(w), tree_root, 0);
            }
        else {// åñëè óçëû - ìåòîäîì choose_node_and_insert íà îïðåäåëåííóþ ãëóáèíó
            for (TreePart* w : forced_reinserted_nodes) {
                choose_node_and_insert(static_cast<Node*>(w), tree_root, deep, 0);
            }
        }
    }

    SplitParameters choose_split_axis_and_index(Node* node) {
        int min_margin(std::numeric_limits<int>::max());
        int distribution_count(max_child_items - 2 * min_child_items + 2);
        SplitParameters params;
        BoundingBox b1, b2;
        for (int axis(0); axis < dimensions; axis++) {//âûáèðàåòñÿ îñü íàèëó÷øåãî ðàñïðåäåëåíèÿ
            for (int i(0); i < 2; i++) {// ïî áîëüøåé èëè ïî ìåíüøåé ãðàíèöå
                axis_type type;
                if (i == 0) {
                    type = axis_type::lower;
                }
                else {
                    type = axis_type::upper;
                }
                std::sort(node->items.begin(), node->items.end(),
                    [&axis, &type](auto& lhs, auto& rhs) {
                    return lhs->box.value_of_axis(axis, type)
                        < rhs->box.value_of_axis(axis, type);
                });
                for (int k(0); k < distribution_count; k++) {// âûáèðàïòñÿ èíäåêñ ëó÷øåãî ðàñïðåäåëåíèÿ
                    int area(0);
                    b1.reset();
                    b2.reset();
                    for (int i = 0; i < min_child_items + k; i++) {
                        b1.stretch(node->items[i]->box);
                    }
                    for (int i = min_child_items + k; i < max_child_items + 1; i++) {
                        b2.stretch(node->items[i]->box);
                    }
                    int margin = b1.margin() + b2.margin();
                    if (margin < min_margin) { // îïðåäåëÿåòñÿ ïî ìèíèìàëüíîé ñóììå ïåðèìåòðîâ
                        min_margin = margin;
                        params.index = k;
                        params.axis = axis;
                        params.type = type;
                    }
                }
            }
        }
        return params;
    }

    void delete_tree(Node* node) {
        if (node->hasleaves) {
            for (size_t i = 0; i < node->items.size(); i++) {
                delete node->items[i];
            }
        }
        else {
            for (size_t i = 0; i < node->items.size(); i++) {
                delete_tree(static_cast<Node*>(node->items[i]));
            }
        }
        delete node;
    }

private:

    std::unordered_set<int> used_deeps;//<используемы границы во время текущей вставки
    Node* tree_root{ nullptr };//<корень дерева
    std::size_t size_{ 0 };//<количество листьев

};