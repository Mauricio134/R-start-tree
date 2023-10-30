#pragma once
#include<vector>
#include<algorithm>
enum class axis_type {
    upper, lower
};

 /*!
 * \brief Шаблонный класс ограничивающей области
 * \param[in] dimensions - количество измерений
 *  */
template<std::size_t dimensions>
struct RStarBoundingBox {
    /*!
 	* \brief Создается экземпляр и заполняется крайними значениями
 	*  */
     RStarBoundingBox() 
         : max_edges(dimensions)
         , min_edges(dimensions) {
         for (size_t axis = 0; axis < dimensions; axis++) {
             max_edges[axis] = std::numeric_limits<int>::min();
             min_edges[axis] = std::numeric_limits<int>::max();
         }
     }
     ~RStarBoundingBox() = default;
     bool operator<(const RStarBoundingBox& rhs) const{
         if (min_edges == rhs.min_edges)
             return max_edges < rhs.max_edges;
         else
             return min_edges < rhs.min_edges;
     }
     bool operator==(const RStarBoundingBox& rhs) const {
         return max_edges == max_edges && min_edges == min_edges;
     }
     bool operator!=(const RStarBoundingBox& rhs)const {
         return !operator==(rhs);
     }

    /*!
 	* \brief Приводит границы области к значениям по умолчанию
 	*  */
     void reset() {
         for (size_t axis = 0; axis < dimensions; axis++) {
             max_edges[axis] = std::numeric_limits<int>::min();
             min_edges[axis] = std::numeric_limits<int>::max();
         }
     }

    /*!
 	* \brief Подбирает границы так, чтобы другая область полностью помещалась в неё
 	* \param[in] other_box - другая область
 	*  */
     void stretch(const RStarBoundingBox<dimensions>& other_box) {
         for (size_t axis = 0; axis < dimensions; axis++) {
             max_edges[axis] = std::max(max_edges[axis], other_box.max_edges[axis]);//Ïîäáèðàåò ãðàíèöû òàê,
             min_edges[axis] = std::min(min_edges[axis], other_box.min_edges[axis]);//÷òîáû äðóãàÿ îáëàñòü ìîãëà ïîìåñòèòüñÿ
         }
     }

    /*!
 	* \brief Возвращает true, если области пересекаются
 	* \param[in] other_box - вторая область
 	*  */
     bool is_intersected(const RStarBoundingBox<dimensions>& other_box) const {
         if (overlap(other_box) > 0)return true;//Âîçâðàùàåò true, åñëè ïåðåñå÷åíèå áîëüøå 0
         else return false;
     }

    /*!
 	* \brief Возвращает полупериметр ограничивающей области
 	*  */   
     int margin() const {
         int ans = 0;
         for (size_t axis = 0; axis < dimensions; axis++) {//Ñ÷èòàåò cóììó ãðàíèö
             ans += max_edges[axis] - min_edges[axis];
         }
         return ans;
     }

    /*!
 	* \brief Возвращает площадь ограничивающей области
 	*  */
     int area() const {//Ñ÷èòàåò ïëîùàäü ïîâåðõíîñòè
         int ans = 1;
         for (size_t axis = 0; axis < dimensions; axis++) {
             ans *= max_edges[axis] - min_edges[axis];
         }
         return ans;
     }

	/*!
 	* \brief Возвращает пересечение двух областей
 	*  */
     int overlap(const RStarBoundingBox<dimensions>& other_box) const {
         int ans = 1;
         for (size_t axis = 0; axis < dimensions; axis++) {
             int x1 = min_edges[axis];// ìåíüøàÿ ãðàíèöà ïåðâîãî
             int x2 = max_edges[axis];// áîëüøàÿ ãðàíèöà ïåðâîãî
             int y1 = other_box.min_edges[axis];//ìåíüøàÿ ãðàíèöà âòîðîãî
             int y2 = other_box.max_edges[axis];//áîëüøàÿ ãðàíèöà âòîðîãî
             
             if (x1 < y1) {
                 if (y1 < x2) {
                     if (y2 < x2) {
                         ans *= (y2 - y1);
                     }
                     else {
                         ans *= (x2 - y1);
                     }
                 }
                 else {
                     return 0;
                 }
             }
             else {
                 if (y2 > x1) {
                     if (y2 > x2) {
                         ans *= (x2 - x1);
                     }
                     else {
                         ans *= (y2 - x1);
                     }
                 }
                 else {
                     return 0; //åñëè íå ïîäîøëî íè â êàêîå óñëîâèå
                 }
                 
             }
         }
         return ans;
     }

     /*!
 	* \brief Возвращает расстояние между центрами двух областей
 	*  */
     double dist_between_centers(const RStarBoundingBox<dimensions>& other_box) const {
         //Ðåçóëüòàò - êâàäðàò ðàññòîÿíèÿ
         //÷òîáû íå òåðÿòü òî÷íîñòü îò êâàäðàòíîãî êîðíÿ
         int ans = 0;
         for (size_t axis = 0; axis < dimensions; axis++) {
             int d = ((max_edges[axis] + min_edges[axis])
                 - (other_box.max_edges[axis] + other_box.min_edges[axis])) / 2;
             ans += d * d;
         }
         return ans;
     }

     int value_of_axis(const int axis, const axis_type type) const {
         if (type == axis_type::lower) {
             return min_edges[axis];
         }
         else {
             return max_edges[axis];
         }
     }
     std::vector<int> max_edges, min_edges; //<границы
};