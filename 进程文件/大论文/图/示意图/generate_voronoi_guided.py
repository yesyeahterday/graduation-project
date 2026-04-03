import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Polygon, FancyArrowPatch
from matplotlib.collections import LineCollection
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import Point, Polygon as ShapelyPolygon, LineString
from shapely.ops import unary_union
import math

def is_line_of_sight(p1, p2, buildings):
    """检查两点之间是否有视距（不被建筑物遮挡）"""
    line = LineString([p1, p2])
    for building in buildings:
        if line.intersects(building):
            return False
    return True

def generate_sensing_area(center, radius, buildings, grid_resolution=50):
    """生成感知区域，考虑建筑物遮挡"""
    x_min, y_min = center[0] - radius, center[1] - radius
    x_max, y_max = center[0] + radius, center[1] + radius
    
    # 生成网格点
    x = np.linspace(x_min, x_max, grid_resolution)
    y = np.linspace(y_min, y_max, grid_resolution)
    X, Y = np.meshgrid(x, y)
    
    points = []
    for i in range(len(x)):
        for j in range(len(y)):
            p = (X[i, j], Y[i, j])
            dist = np.sqrt((p[0] - center[0])**2 + (p[1] - center[1])**2)
            if dist <= radius:
                if is_line_of_sight(center, p, buildings):
                    points.append(p)
    
    return points

def calculate_polygon_centroid(vertices):
    """计算多边形的质心"""
    if len(vertices) < 3:
        return None
    
    x = [v[0] for v in vertices]
    y = [v[1] for v in vertices]
    
    # 使用Shoelace公式计算质心
    area = 0
    cx = 0
    cy = 0
    
    for i in range(len(vertices)):
        j = (i + 1) % len(vertices)
        cross = x[i] * y[j] - x[j] * y[i]
        area += cross
        cx += (x[i] + x[j]) * cross
        cy += (y[i] + y[j]) * cross
    
    if abs(area) < 1e-10:
        return None
    
    area /= 2
    cx /= (6 * area)
    cy /= (6 * area)
    
    return (cx, cy)

def find_neighbors(selected_idx, drones, radius, buildings):
    """找到距离r范围内的邻居节点（考虑视距）"""
    selected = drones[selected_idx]
    neighbors = []
    
    for i, drone in enumerate(drones):
        if i == selected_idx:
            continue
        
        dist = np.sqrt((drone[0] - selected[0])**2 + (drone[1] - selected[1])**2)
        if dist <= radius:
            if is_line_of_sight(selected, drone, buildings):
                neighbors.append(i)
    
    return neighbors

def draw_perpendicular_bisector(p1, p2, bounds, ax):
    """绘制两点之间的垂直平分线"""
    # 计算中点
    mid_x = (p1[0] + p2[0]) / 2
    mid_y = (p1[1] + p2[1]) / 2
    
    # 计算方向向量
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    
    # 如果两点重合，不绘制
    if abs(dx) < 1e-10 and abs(dy) < 1e-10:
        return
    
    # 垂直方向（旋转90度）
    perp_dx = -dy
    perp_dy = dx
    
    # 归一化
    length = np.sqrt(perp_dx**2 + perp_dy**2)
    if length < 1e-10:
        return
    perp_dx /= length
    perp_dy /= length
    
    # 计算线段长度（延伸到边界外）
    x_min, y_min, x_max, y_max = bounds
    max_dist = max(x_max - x_min, y_max - y_min) * 2
    
    # 计算线段端点（在边界外）
    start_x = mid_x - perp_dx * max_dist
    start_y = mid_y - perp_dy * max_dist
    end_x = mid_x + perp_dx * max_dist
    end_y = mid_y + perp_dy * max_dist
    
    # 使用Liang-Barsky算法裁剪到边界
    def clip_line(x0, y0, x1, y1, x_min, y_min, x_max, y_max):
        """裁剪线段到矩形边界内"""
        dx = x1 - x0
        dy = y1 - y0
        
        if abs(dx) < 1e-10 and abs(dy) < 1e-10:
            if x_min <= x0 <= x_max and y_min <= y0 <= y_max:
                return [(x0, y0), (x1, y1)]
            return None
        
        t0, t1 = 0.0, 1.0
        
        # 检查各个边界
        for edge in range(4):
            if edge == 0:   # 左边界
                p, q = -dx, x0 - x_min
            elif edge == 1: # 右边界
                p, q = dx, x_max - x0
            elif edge == 2: # 下边界
                p, q = -dy, y0 - y_min
            else:           # 上边界
                p, q = dy, y_max - y0
            
            if abs(p) < 1e-10:
                if q < 0:
                    return None
            else:
                r = q / p
                if p < 0:
                    if r > t1:
                        return None
                    elif r > t0:
                        t0 = r
                else:
                    if r < t0:
                        return None
                    elif r < t1:
                        t1 = r
        
        # 计算裁剪后的端点
        new_x0 = x0 + t0 * dx
        new_y0 = y0 + t0 * dy
        new_x1 = x0 + t1 * dx
        new_y1 = y0 + t1 * dy
        
        return [(new_x0, new_y0), (new_x1, new_y1)]
    
    clipped = clip_line(start_x, start_y, end_x, end_y, x_min, y_min, x_max, y_max)
    
    if clipped:
        # 绘制垂直平分线（黑色虚线）
        ax.plot(
            [clipped[0][0], clipped[1][0]],
            [clipped[0][1], clipped[1][1]],
            linestyle="--",
            color="black",
            linewidth=1.5,
            alpha=0.8,
            zorder=5,
        )

def get_voronoi_cell(selected_idx, neighbors, drones, bounds):
    """获取选中节点在Voronoi图中的晶胞"""
    # 构建包含选中节点和邻居节点的点集
    points = [drones[selected_idx]]
    for idx in neighbors:
        points.append(drones[idx])
    
    # 如果没有邻居，返回整个区域
    if len(points) == 1:
        x_min, y_min, x_max, y_max = bounds
        return [[x_min, y_min], [x_max, y_min], [x_max, y_max], [x_min, y_max]]
    
    # 添加边界点以避免无限区域
    x_min, y_min, x_max, y_max = bounds
    margin = max(x_max - x_min, y_max - y_min) * 0.5
    boundary_points = [
        [x_min - margin, y_min - margin],
        [x_max + margin, y_min - margin],
        [x_max + margin, y_max + margin],
        [x_min - margin, y_max + margin]
    ]
    all_points = boundary_points + points
    
    # 计算Voronoi图
    try:
        vor = Voronoi(all_points)
        
        # 找到选中节点对应的区域（索引为len(boundary_points)，即4）
        region_idx = len(boundary_points)
        if region_idx >= len(vor.point_region):
            return None
        
        region_id = vor.point_region[region_idx]
        if region_id == -1 or region_id >= len(vor.regions):
            return None
        
        region = vor.regions[region_id]
        
        # 获取区域顶点
        vertices = []
        for vertex_idx in region:
            if vertex_idx != -1 and vertex_idx < len(vor.vertices):
                vertices.append(tuple(vor.vertices[vertex_idx]))
        
        # 裁剪到边界内
        if len(vertices) >= 3:
            try:
                polygon = ShapelyPolygon(vertices)
                bounds_poly = ShapelyPolygon([
                    [x_min, y_min], [x_max, y_min],
                    [x_max, y_max], [x_min, y_max]
                ])
                clipped = polygon.intersection(bounds_poly)
                if clipped.geom_type == 'Polygon' and not clipped.is_empty:
                    vertices = list(clipped.exterior.coords[:-1])
                elif clipped.geom_type == 'MultiPolygon':
                    # 取最大的多边形
                    largest = max(clipped.geoms, key=lambda p: p.area)
                    vertices = list(largest.exterior.coords[:-1])
            except:
                pass
        
        return vertices if len(vertices) >= 3 else None
    except:
        return None

def main():
    # 参数设置
    area_size = 50  # 正方形区域边长
    num_buildings = 6  # 建筑物数量
    num_drones = 7  # 无人机数量
    sensing_radius = 10  # 感知半径（相应缩小）
    neighbor_radius = 30  # 邻居节点选择半径（相应缩小）
    selected_drone_idx = 1  # 选中的无人机索引
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_xlim(0, area_size)
    ax.set_ylim(0, area_size)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # 1. 绘制正方形区域（黑边1.5磅）
    square = Rectangle((0, 0), area_size, area_size, 
                      linewidth=1.5, edgecolor='black', 
                      facecolor='white', fill=True, zorder=0)
    ax.add_patch(square)
    
    # 2. 生成随机建筑物
    np.random.seed(42)  # 设置随机种子以便复现
    buildings = []
    building_patches = []
    min_distance = 5  # 建筑物之间的最小距离
    
    def check_building_distance(new_x, new_y, new_width, new_height, existing_buildings, min_dist):
        """检查新建筑物与已有建筑物之间的距离"""
        new_center_x = new_x + new_width / 2
        new_center_y = new_y + new_height / 2
        
        for existing in existing_buildings:
            # 获取已有建筑物的中心
            existing_bounds = existing.bounds
            existing_center_x = (existing_bounds[0] + existing_bounds[2]) / 2
            existing_center_y = (existing_bounds[1] + existing_bounds[3]) / 2
            
            # 计算中心距离
            dist = np.sqrt((new_center_x - existing_center_x)**2 + 
                          (new_center_y - existing_center_y)**2)
            
            # 计算最小允许距离（考虑建筑物大小）
            existing_width = existing_bounds[2] - existing_bounds[0]
            existing_height = existing_bounds[3] - existing_bounds[1]
            min_allowed_dist = (new_width + new_height + existing_width + existing_height) / 4 + min_dist
            
            if dist < min_allowed_dist:
                return False
        return True
    
    for i in range(num_buildings):
        max_attempts = 100
        attempt = 0
        placed = False
        
        while not placed and attempt < max_attempts:
            # 随机生成建筑物位置和大小，尺寸差异更大
            # 小建筑物：2-5，大建筑物：8-15
            if np.random.random() < 0.5:
                # 小建筑物
                width = np.random.uniform(2, 5)
                height = np.random.uniform(2, 5)
            else:
                # 大建筑物
                width = np.random.uniform(8, 15)
                height = np.random.uniform(8, 15)
            
            x = np.random.uniform(2, area_size - width - 2)
            y = np.random.uniform(2, area_size - height - 2)
            
            # 确保不超出边界
            if x + width > area_size:
                width = area_size - x - 1
            if y + height > area_size:
                height = area_size - y - 1
            
            if width < 2 or height < 2:
                attempt += 1
                continue
            
            # 检查与已有建筑物的距离
            if i == 0 or check_building_distance(x, y, width, height, buildings, min_distance):
                # 使用指定的黄色 #FFD966
                building = Rectangle(
                    (x, y),
                    width,
                    height,
                    facecolor="#FFD966",
                    edgecolor="black",
                    linewidth=2,
                    alpha=0.85,
                )
                ax.add_patch(building)
                building_patches.append(building)
                
                # 创建Shapely多边形用于视距检测
                building_poly = ShapelyPolygon([
                    [x, y], [x + width, y],
                    [x + width, y + height], [x, y + height]
                ])
                buildings.append(building_poly)
                placed = True
            else:
                attempt += 1

    # 将从左到右排序后的第二个建筑物移动到指定位置附近 (10, 30)
    if len(building_patches) >= 2:
        # 按 x 坐标从小到大排序，取第二个
        sorted_indices = sorted(
            range(len(building_patches)),
            key=lambda idx: building_patches[idx].get_x(),
        )
        second_idx = sorted_indices[1]
        patch = building_patches[second_idx]

        width = patch.get_width()
        height = patch.get_height()

        # 目标坐标（略作裁剪，确保在区域内）
        target_x = 10.0
        target_y = 30.0
        new_x = max(1.0, min(target_x, area_size - width - 1.0))
        new_y = max(1.0, min(target_y, area_size - height - 1.0))

        patch.set_x(new_x)
        patch.set_y(new_y)

        # 同步更新对应的 Shapely 多边形
        bpoly = buildings[second_idx]
        buildings[second_idx] = ShapelyPolygon(
            [
                [new_x, new_y],
                [new_x + width, new_y],
                [new_x + width, new_y + height],
                [new_x, new_y + height],
            ]
        )
    
    # 3. 生成无人机位置（缩小分布范围）
    drones = []
    for i in range(num_drones):
        # 确保无人机不在建筑物内
        while True:
            x = np.random.uniform(2, area_size - 2)
            y = np.random.uniform(2, area_size - 2)
            point = Point(x, y)
            
            in_building = False
            for building in buildings:
                if building.contains(point):
                    in_building = True
                    break
            
            if not in_building:
                drones.append([x, y])
                break
    
    # 绘制无人机（黑色质点）
    for i, drone in enumerate(drones):
        if i == selected_drone_idx:
            # 选中的节点用绿色
            ax.plot(drone[0], drone[1], 'o', color='green', 
                   markersize=10, markeredgecolor='black', 
                   markeredgewidth=1.5, label='选中节点', zorder=10)
        else:
            ax.plot(drone[0], drone[1], 'o', color='black', 
                   markersize=8, zorder=10)
    
    # 4. 绘制感知区域（更明显的蓝色，考虑非视距）
    for i, drone in enumerate(drones):
        sensing_points = generate_sensing_area(drone, sensing_radius, buildings)
        if sensing_points:
            x_coords = [p[0] for p in sensing_points]
            y_coords = [p[1] for p in sensing_points]
            ax.scatter(
                x_coords,
                y_coords,
                c="#66B3FF",        # 更饱和的蓝色
                s=6,               # 增大点的尺寸
                alpha=0.7,         # 提高不透明度
                edgecolors="none",
                zorder=1,
            )
    
    # 5. Voronoi算法示意
    # 找到邻居节点
    neighbors = find_neighbors(selected_drone_idx, drones, neighbor_radius, buildings)
    
    # 绘制邻居节点（蓝色）
    if neighbors:
        ax.plot([drones[idx][0] for idx in neighbors], 
               [drones[idx][1] for idx in neighbors], 'o', 
               color='blue', markersize=8, 
               markeredgecolor='black', markeredgewidth=1,
               label='邻居节点', zorder=10)
        
        # 绘制与每个邻居节点的垂直平分线
        selected_drone = drones[selected_drone_idx]
        for neighbor_idx in neighbors:
            neighbor_drone = drones[neighbor_idx]
            draw_perpendicular_bisector(selected_drone, neighbor_drone, 
                                       (0, 0, area_size, area_size), ax)
    
    # 获取Voronoi晶胞
    bounds = (0, 0, area_size, area_size)
    cell_vertices = get_voronoi_cell(selected_drone_idx, neighbors, drones, bounds)
    
    if cell_vertices and len(cell_vertices) > 2:
        # 原始晶胞多边形
        cell_poly = ShapelyPolygon(cell_vertices)

        # 先根据建筑物从晶胞中扣除建筑体本身，作为“可见区域”的基础
        occluded_polys = []
        for bpoly in buildings:
            if cell_poly.intersects(bpoly):
                inter = cell_poly.intersection(bpoly)
                if not inter.is_empty:
                    occluded_polys.append(inter)

        if occluded_polys:
            occluded_union = unary_union(occluded_polys)
            visible_poly = cell_poly.difference(occluded_union)
        else:
            visible_poly = cell_poly

        # 绘制原始晶胞边界（黑色实线）
        exterior = cell_poly.exterior.coords
        xs = [p[0] for p in exterior]
        ys = [p[1] for p in exterior]
        ax.plot(xs, ys, color="black", linewidth=2, linestyle="-", alpha=0.9, zorder=7)

        # 基于建筑物生成“遮挡扇区”的两根虚线 / 实线射线
        selected_drone = drones[selected_drone_idx]
        sx, sy = selected_drone
        x_min, y_min, x_max, y_max = bounds
        max_dist = max(x_max - x_min, y_max - y_min) * 2

        for bpoly in buildings:
            if not cell_poly.intersects(bpoly):
                continue

            # 取建筑物的外轮廓顶点，按与选中节点的角度排序
            vertices = list(bpoly.exterior.coords)[:-1]
            angles = []
            for vx, vy in vertices:
                dx = vx - sx
                dy = vy - sy
                if abs(dx) < 1e-8 and abs(dy) < 1e-8:
                    continue
                angle = math.atan2(dy, dx)
                angles.append((angle, (vx, vy)))

            if len(angles) < 2:
                continue

            angles.sort(key=lambda t: t[0])
            # 选择最小和最大的两个角度，近似表示遮挡扇区边界
            extreme = [angles[0], angles[-1]]

            for ang, _ in extreme:
                dir_x = math.cos(ang)
                dir_y = math.sin(ang)
                # 从节点向外发射射线
                ray_end = (sx + dir_x * max_dist, sy + dir_y * max_dist)
                ray = LineString([(sx, sy), ray_end])

                # 与晶胞的交点（确定射线在晶胞内的最远端）
                inter_cell = ray.intersection(cell_poly)
                if inter_cell.is_empty:
                    continue
                cell_pts = []
                if inter_cell.geom_type == "Point":
                    cell_pts = [inter_cell]
                elif inter_cell.geom_type in ("MultiPoint", "GeometryCollection"):
                    cell_pts = [g for g in inter_cell.geoms if g.geom_type == "Point"]
                elif inter_cell.geom_type == "LineString":
                    cell_pts = [Point(c) for c in inter_cell.coords]
                if not cell_pts:
                    continue
                far_pt = max(
                    cell_pts,
                    key=lambda p: (p.x - sx) ** 2 + (p.y - sy) ** 2,
                )

                # 与建筑物的交点（确定遮挡开始位置）
                inter_build = ray.intersection(bpoly)
                build_pts = []
                if not inter_build.is_empty:
                    if inter_build.geom_type == "Point":
                        build_pts = [inter_build]
                    elif inter_build.geom_type in ("MultiPoint", "GeometryCollection"):
                        build_pts = [g for g in inter_build.geoms if g.geom_type == "Point"]
                    elif inter_build.geom_type == "LineString":
                        build_pts = [Point(c) for c in inter_build.coords]

                if build_pts:
                    near_build = min(
                        build_pts,
                        key=lambda p: (p.x - sx) ** 2 + (p.y - sy) ** 2,
                    )
                    d_cell = (far_pt.x - sx) ** 2 + (far_pt.y - sy) ** 2
                    d_build = (near_build.x - sx) ** 2 + (near_build.y - sy) ** 2

                    # 从节点到建筑边界：视距区域，用虚线
                    ax.plot(
                        [sx, near_build.x],
                        [sy, near_build.y],
                        linestyle="--",
                        color="black",
                        linewidth=2.0,
                        alpha=0.7,
                        zorder=6,
                    )

                    # 从建筑边界到晶胞边界：被建筑遮挡的背面，用实线
                    if d_build < d_cell:
                        ax.plot(
                            [near_build.x, far_pt.x],
                            [near_build.y, far_pt.y],
                            linestyle="-",
                            color="black",
                            linewidth=2.0,
                            alpha=0.7,
                            zorder=6,
                        )
                else:
                    # 与建筑不发生有效相交时，全程视距，用实线
                    ax.plot(
                        [sx, far_pt.x],
                        [sy, far_pt.y],
                        linestyle="-",
                        color="black",
                        linewidth=2.0,
                        alpha=0.7,
                        zorder=6,
                    )

        # 计算晶胞质心（以可见区域为主，若为空则退回原始晶胞）
        centroid_source = visible_poly if (not visible_poly.is_empty) else cell_poly
        if centroid_source.is_empty:
            centroid_point = Point(
                sum(p[0] for p in cell_vertices) / len(cell_vertices),
                sum(p[1] for p in cell_vertices) / len(cell_vertices),
            )
        else:
            if isinstance(centroid_source, ShapelyPolygon):
                centroid_point = centroid_source.centroid
            else:
                # MultiPolygon：选择面积最大的那一个
                largest = max(centroid_source.geoms, key=lambda g: g.area)
                centroid_point = largest.centroid

        centroid = (centroid_point.x, centroid_point.y)

        # 绘制质心（三角形标记）
        ax.plot(
            centroid[0],
            centroid[1],
            "^",
            color="red",
            markersize=12,
            markeredgecolor="black",
            markeredgewidth=1.5,
            label="质心",
            zorder=10,
        )

        # 绘制从节点指向质心的箭头
        dx = centroid[0] - selected_drone[0]
        dy = centroid[1] - selected_drone[1]
        dist = np.sqrt(dx**2 + dy**2)
        if dist > 0:
            # 箭头起点稍微远离节点，终点稍微远离质心
            start_offset = 1.5
            end_offset = 1.5
            start_x = selected_drone[0] + dx * start_offset / dist
            start_y = selected_drone[1] + dy * start_offset / dist
            end_x = centroid[0] - dx * end_offset / dist
            end_y = centroid[1] - dy * end_offset / dist

            arrow = FancyArrowPatch(
                (start_x, start_y),
                (end_x, end_y),
                arrowstyle="->",
                mutation_scale=20,
                linewidth=2,
                color="red",
                zorder=10,
            )
            ax.add_patch(arrow)
    
    # 添加图例
    ax.legend(loc='upper right', fontsize=9)
    
    # 保存为SVG
    plt.tight_layout()
    plt.savefig('voronoi示意图.svg', format='svg', dpi=300, bbox_inches='tight')
    print("示意图已保存为 voronoi示意图.svg")
    
    plt.show()

if __name__ == "__main__":
    main()